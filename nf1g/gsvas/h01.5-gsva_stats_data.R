source("libs.R")
library(dtplyr)
library(purrr)
library(tidyverse)
library(crayon)



hgs0<-readRDS("timex/ds/hu-msig_all-gsva-values.rds")

sample_info<-readRDS("nf1g/ds/v10-per_sample_updated.rds")





patho_cat_names<-sample_info%>%
  select(patho_cat_name)%>%
  unique()%>%
  pull()%>%
  as.character()


patho_cat_names[1:5]


for (p in 1:length(patho_cat_names)){
# p<-3


file_name_ext<-patho_cat_names[p]
  
  
patho_cat_names[p]



samples_subset_info<-sample_info%>%
  filter(patho_cat_name==patho_cat_names[p])

group_factors<-samples_subset_info%>%
  select(resultant_geno)%>%unique()%>%
  pull(resultant_geno)


group_a_samples<-samples_subset_info%>%
  filter(resultant_geno==group_factors[1])%>%
  pull(mouse_num)

group_b_samples<-samples_subset_info%>%
  filter(resultant_geno==group_factors[2])%>%
  pull(mouse_num)



if(length(group_a_samples)>=2&length(group_b_samples)>=2){



sub_matrix<-hgs0[,c(group_a_samples, group_b_samples)]



# Ensure sample names match matrix column names
valid_a_samples <- intersect(group_a_samples, colnames(sub_matrix))
valid_b_samples <- intersect(group_b_samples, colnames(sub_matrix))

# Extract the relevant subsets of the matrix for each group
group_a_matrix <- sub_matrix[, valid_a_samples, drop = FALSE]
group_b_matrix <- sub_matrix[, valid_b_samples, drop = FALSE]



# Initialize vectors to store results
fold_changes <- numeric(nrow(sub_matrix))
p_values_ttest <- numeric(nrow(sub_matrix))
p_values_lm <- numeric(nrow(sub_matrix))
mean_total<-numeric(nrow(sub_matrix))
mean_a <- numeric(nrow(sub_matrix))
mean_b <- numeric(nrow(sub_matrix))
std_a <- numeric(nrow(sub_matrix))
std_b <- numeric(nrow(sub_matrix))
analysis_method <- character(nrow(sub_matrix))  # Store method used for p-value




# Apply t-test and linear model for each pathway (row)
for (i in 1:nrow(sub_matrix)) {
  # Get z-scores for this pathway across samples
  a_values <- group_a_matrix[i, , drop = FALSE]
  b_values <- group_b_matrix[i, , drop = FALSE]
  
  # Combine values into a single vector for fitting the linear model
  combined_values <- c(a_values, b_values)
  group_labels <- rep(c(group_factors[1], group_factors[2]), times = c(length(valid_a_samples), length(valid_b_samples)))
  
  # Prepare data for linear model
  lm_data <- data.frame(z_score = combined_values, Group = factor(group_labels))
  
  # Perform t-test for this pathway
  ttest_result <- t.test(a_values, b_values)
  p_values_ttest[i] <- ttest_result$p.value
  
  # Store analysis method for t-test
  analysis_method[i] <- "t-test"
  
  # Fit the linear model: z_score ~ Group (group comparison)
  lm_model <- lm(z_score ~ Group, data = lm_data)
  
  # Store results for linear model
  p_values_lm[i] <- summary(lm_model)$coefficients[2, 4]  # p-value for group comparison
  fold_changes[i] <- coef(lm_model)[2]  # Mean difference (group B vs group A)
  
  # Store analysis method for linear model if p-value is more significant
  if (p_values_lm[i] < p_values_ttest[i]) {
    analysis_method[i] <- "lm"
  }
  
  # Compute mean for each group
  mean_a[i] <- mean(a_values, na.rm = TRUE)
  mean_b[i] <- mean(b_values, na.rm = TRUE)
  
  
  mean_total[i] <- mean(c(a_values,b_values), na.rm = TRUE)
  
  
  # Compute standard deviation for each group
  std_a[i] <- sd(a_values, na.rm = TRUE)
  std_b[i] <- sd(b_values, na.rm = TRUE)
  
  
  if (i %% 10 == 0) {
    cat(blue(i), "out of", green(nrow(sub_matrix)), " ", red(round(i/nrow(sub_matrix)*100, 2)), "%", "\n")
  }
  
  
}





# Create results tibble
results <- tibble(
  Pathway = rownames(sub_matrix),
  p_value_ttest = p_values_ttest,
  p_value_lm = p_values_lm,
  fold_change = fold_changes,
  mean_total=mean_total,
  mean_group_a = mean_a,
  mean_group_b = mean_b,
  std_group_a = std_a,
  std_group_b = std_b,
  analysis_method = analysis_method,
  diff = mean_b-mean_a,
  abs_diff=abs(diff)
) %>%
  arrange(p_value_lm)  # Sort by linear model p-value for significance





# Create the output list with metadata
analysis_output <- list(
  results = results,
  metadata = list(
    date = Sys.time(),
    group_comparison = paste(group_factors[1], "vs", group_factors[2]),
    group_a = paste(group_factors[1]),
    group_b = paste(group_factors[2]),
    group_a_samples = valid_a_samples,
    group_b_samples = valid_b_samples,
    src=paste0("src: ",
               
               rstudioapi::getSourceEditorContext()$path%>%
                 sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                 sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
               
               " at ", 
               
               lubridate::round_date(Sys.time(), "second")
    )
  )
)





analysis_output$results2 <- analysis_output$results %>%
  mutate(min_pval = if_else(p_value_lm<=p_value_ttest, p_value_lm, p_value_ttest))




saveRDS(analysis_output, paste0("nf1g/gsvas/ds/hu-gsva_pathway_stats_",
        
        file_name_ext, 
        
        ".rds")

)








}else{
  

  
  

  sub_matrix<-hgs0[,c(group_a_samples, group_b_samples)]
  
  
  
  # Ensure sample names match matrix column names
  valid_a_samples <- intersect(group_a_samples, colnames(sub_matrix))
  valid_b_samples <- intersect(group_b_samples, colnames(sub_matrix))
  
  # Extract the relevant subsets of the matrix for each group
  group_a_matrix <- sub_matrix[, valid_a_samples, drop = FALSE]
  group_b_matrix <- sub_matrix[, valid_b_samples, drop = FALSE]
  
  
  
  # Initialize vectors to store results
  fold_changes <- numeric(nrow(sub_matrix))
  #p_values_ttest <- numeric(nrow(sub_matrix))
  #p_values_lm <- numeric(nrow(sub_matrix))
  mean_total<-numeric(nrow(sub_matrix))
  mean_a <- numeric(nrow(sub_matrix))
  mean_b <- numeric(nrow(sub_matrix))
  
  #std_a <- numeric(nrow(sub_matrix))
  #std_b <- numeric(nrow(sub_matrix))
  
  analysis_method <- character(nrow(sub_matrix))  # Store method used for p-value
  
  
  
  
  # Apply t-test and linear model for each pathway (row)
  for (i in 1:nrow(sub_matrix)) {
    # Get z-scores for this pathway across samples
    a_values <- group_a_matrix[i, , drop = FALSE]
    b_values <- group_b_matrix[i, , drop = FALSE]
    
    # Combine values into a single vector for fitting the linear model
    combined_values <- c(a_values, b_values)
    group_labels <- rep(c(group_factors[1], group_factors[2]), times = c(length(valid_a_samples), length(valid_b_samples)))
    
    # Prepare data for linear model
    #lm_data <- data.frame(z_score = combined_values, Group = factor(group_labels))
    
    # Perform t-test for this pathway
    #ttest_result <- t.test(a_values, b_values)
    #p_values_ttest[i] <- ttest_result$p.value
    
    # Store analysis method for t-test
    #analysis_method[i] <- "t-test"
    
    # Fit the linear model: z_score ~ Group (group comparison)
    #lm_model <- lm(z_score ~ Group, data = lm_data)
    
    # Store results for linear model
    #p_values_lm[i] <- summary(lm_model)$coefficients[2, 4]  # p-value for group comparison
    #fold_changes[i] <- coef(lm_model)[2]  # Mean difference (group B vs group A)
    
    # Store analysis method for linear model if p-value is more significant
    #if (p_values_lm[i] < p_values_ttest[i]) {
     # analysis_method[i] <- "lm"
    #}
    
    # Compute mean for each group
    mean_a[i] <- mean(a_values, na.rm = TRUE)
    mean_b[i] <- mean(b_values, na.rm = TRUE)
    
    
    mean_total[i] <- mean(c(a_values,b_values), na.rm = TRUE)
    
    
    # Compute standard deviation for each group
    std_a[i] <- sd(a_values, na.rm = TRUE)
    std_b[i] <- sd(b_values, na.rm = TRUE)
    
    
    if (i %% 10 == 0) {
      cat(blue(i), "out of", green(nrow(sub_matrix)), " ", red(round(i/nrow(sub_matrix)*100, 2)), "%", "\n")
    }
    
    
  }
  
  
  
  
  
  # Create results tibble
  results <- tibble(
    Pathway = rownames(sub_matrix),
    #p_value_ttest = p_values_ttest,
    #p_value_lm = p_values_lm,
    #fold_change = fold_changes,
    mean_total=mean_total,
    mean_group_a = mean_a,
    mean_group_b = mean_b,
    std_group_a = std_a,
    std_group_b = std_b,
    #analysis_method = analysis_method,
    diff = mean_b-mean_a,
    abs_diff=abs(diff)
  )
    #arrange(p_value_lm)  # Sort by linear model p-value for significance
  
  
  
  
  
  # Create the output list with metadata
  analysis_output <- list(
    results = results,
    metadata = list(
      date = Sys.time(),
      group_comparison = paste(group_factors[1], "vs", group_factors[2]),
      group_a = paste(group_factors[1]),
      group_b = paste(group_factors[2]),
      group_a_samples = valid_a_samples,
      group_b_samples = valid_b_samples,
      src=paste0("src: ",
                 
                 rstudioapi::getSourceEditorContext()$path%>%
                   sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                   sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                 
                 " at ", 
                 
                 lubridate::round_date(Sys.time(), "second")
      )
    )
  )
  
  
  
  
  
  analysis_output$results2 <- analysis_output$results 
    #mutate(min_pval = if_else(p_value_lm<=p_value_ttest, p_value_lm, p_value_ttest))
  
  
  
  
  saveRDS(analysis_output, paste0("nf1g/gsvas/ds/hu-gsva_pathway_stats_a86_",
                                  
                                  file_name_ext, 
                                  
                                  ".rds")
          
  )
  

  
}




}





