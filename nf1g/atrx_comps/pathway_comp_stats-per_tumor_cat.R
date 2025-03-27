source("libs.R")
library(tidyverse)
library(dtplyr)
library(purrr)
standard_error <- function(x) {
  sd(x) / sqrt(length(x))
}

#get in metadata for samples.
all_sample_info<-readRDS("nf1g/ds/v10-per_sample_updated.rds")%>%
  filter(coh==1|coh==2)

#filter for cohorts we're after.
tumor_types<-all_sample_info%>%
  filter(patho_cat!="4")%>%
  select(patho_cat_name)%>%
  unique()%>%
  pull()%>%
  as.character()

tm<-3

tum<-tumor_types[tm]

sample_info<-all_sample_info%>%
  filter(patho_cat_name==tum)

#get simplest naming cohort map to use.
name_map<-sample_info%>%
  select(sample_id, mouse_num, coh, patho_cat_name)

#get mouse nums for each set.
mouse_nums_coh1<-sample_info%>%filter(coh==1)%>%pull(mouse_num)
mouse_nums_coh2<-sample_info%>%filter(coh==2)%>%pull(mouse_num)
mouse_nums<-sample_info%>%pull(mouse_num)


#read in dataset of pathway enrichment scores.
#subset by only mouse_nums we want.
gsva_broad<-readRDS("nf1g/ds/gsva_all_u.rds")

timex_names <- readRDS("nf1g/ds/gsva_sig_u-onco.rds") %>%
  mutate(naked = mapply(function(s, sig) gsub(paste0(s, "_"), "", sig), Set, Signature))%>%
  mutate(fullname=paste0(toupper(Set), "_", naked))%>%
  mutate(fullname = str_replace(fullname, "/", "v."))%>%
  rename(name1=Signature, name2=fullname)%>%
  select(name1,name2)

gsva_timex<-readRDS("nf1g/ds/gsva_u-onco.rds")


rownames(gsva_timex) <- timex_names$name2[match(rownames(gsva_timex), timex_names$name1)]


gsva<-rbind(gsva_broad, gsva_timex)[,mouse_nums]


#convert to a df for other uses.
gsva_df<-gsva%>%as.data.frame()%>%select(all_of(mouse_nums))%>%
  rownames_to_column("pathway")


#pivot gsva values to enable analysis.
pathways_diff <- gsva_df %>%
  pivot_longer(cols = -pathway, names_to = "mouse_num", values_to = "gsea_val")%>%
  left_join(name_map)
  

####generate group comparison p values.
######################## unpaired ############

# Extract unique pathway names
pathway_names <- rownames(gsva)

# Initialize a tibble to store results
all_ps_unpaired <- tibble(pn = character(), pv = numeric())


# Pre-filter the data for coh == 1 and coh == 2 once
coh1_data <- pathways_diff %>% filter(coh == 1) %>% select(pathway, mouse_num, gsea_val)
coh2_data <- pathways_diff %>% filter(coh == 2) %>% select(pathway, mouse_num, gsea_val)

# Convert to matrices (gsea_val only)
coh1_matrix <- coh1_data %>% spread(key = pathway, value = gsea_val) %>% select(-mouse_num) %>% as.matrix()
coh2_matrix <- coh2_data %>% spread(key = pathway, value = gsea_val) %>% select(-mouse_num) %>% as.matrix()

# Initialize a vector to store p-values
p_values <- numeric(length(pathway_names))

pathway_names%>%as_tibble()%>%view()



# Loop through each pathway using matrix indexing
for (i in seq_along(pathway_names)) {
  pn <- pathway_names[i]
  
  # Get the column index for the current pathway
  pathway_index <- which(colnames(coh1_matrix) == pn)
  
  # Check if the pathway exists in both matrices
  if (length(pathway_index) == 0) {
    stop(paste("Pathway", pn, "not found in matrices"))
  }
  
  # Count non-NA values in each cohort
  n1 <- sum(!is.na(coh1_matrix[, pathway_index]))
  n2 <- sum(!is.na(coh2_matrix[, pathway_index]))
  
  # Stop if there are not enough observations
  if (n1 < 2 | n2 < 2) {
    stop(paste("Not enough observations for pathway:", pn, " | Coh1:", n1, " | Coh2:", n2))
  }
  
  # Perform the t-test
  pv <- t.test(coh1_matrix[, pathway_index], coh2_matrix[, pathway_index], paired = FALSE)$p.value
  
  # Store the result
  p_values[i] <- pv
  
  # Print progress every 10 iterations
  remaining <- length(pathway_names) - i
  if (remaining %% 10 == 0) {
    print(paste(remaining, "remaining"))
  }
}

# Create the final tibble with pathway names and p-values
all_ps_unpaired <- tibble(pn = pathway_names, pv = p_values)




#now generate group stats.
###########################################
mouse_nums_list <- list(
  all_mice = mouse_nums, 
  coh1_mice = mouse_nums_coh1, 
  coh2_mice = mouse_nums_coh2
)

# Pre-allocate the list to store summary statistics
stats_list <- list()

# # Loop through each mouse group
for (m in 1:length(mouse_nums_list)) {


  list_name <- names(mouse_nums_list)[[m]]
  
  # Select the mice corresponding to the current group
  mice <- mouse_nums_list[[m]]
  
  # # Convert mice indices to numeric
  # mice <- as.numeric(mice)  # Ensures the indices are numeric
  # 
  
  # Ensure all indices are within bounds of gsva_matrix
 
  
  # Extract the data for the current group, ensuring no out-of-bounds indices
  group_data <- gsva[, mice, drop = FALSE]  # Use column indices
  
  # Compute summary statistics for the group (min, max, mean, sd, se) column-wise
  min_values <- apply(group_data, 1, min)  # Apply on rows (1 indicates rows)
  max_values <- apply(group_data, 1, max)
  mean_values <- apply(group_data, 1, mean)
  sd_values <- apply(group_data, 1, sd)
  se_values <- apply(group_data, 1, function(x) sd(x) / sqrt(length(x)))
  
  # Combine the results into a data frame for the current group
  summary_stats <- cbind(min_value = min_values,
                         max_value = max_values,
                         mean_value = mean_values,
                         sd_value = sd_values,
                         se_value = se_values)
  
  # Store the results in the list with the appropriate name
  stats_list[[list_name]] <- summary_stats
  
  print(paste("Finished group", m))
}


stats_list_df <- lapply(stats_list, function(mat) {
  df <- as.data.frame(mat)%>%
    rownames_to_column("pn")  # Add pathway as a column based on row names

})



stats_tibble2 <- stats_list_df %>%
  map2(names(stats_list_df), function(df, name) {
    df %>%
      rename_with(~ paste0(name, "_", .), -pn)  # Rename all columns except 'pathway'
  }) %>%
  reduce(function(x, y) left_join(x, y, by = "pn"))%>%
  
  mutate(fc=coh2_mice_mean_value/coh1_mice_mean_value)



stats_tibble3<-all_ps_unpaired%>%
  left_join(stats_tibble2)%>%
  relocate(pn, pv, fc)



saveRDS(stats_tibble3, paste0("nf1g/atrx_comps/ds/atrx_comps-gsva", tum, ".rds"))







