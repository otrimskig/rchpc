source("libs.R")
library(dtplyr)
library(purrr)
library(tidyverse)
library(crayon)



hgs0<-readRDS("timex/ds/hu-msig_all-gsva-values.rds")

sample_info<-readRDS("nf1g/ds/v10-per_sample_updated.rds")

samples_subset_info<-sample_info%>%
  filter(resultant_geno!="nf1 wt; pten wt; ink wt; atrx wt")
  


sub_matrix<-hgs0[,samples_subset_info%>%pull(mouse_num)]



# Initialize vectors to store results
mean_total<-numeric(nrow(sub_matrix))
std_total <- numeric(nrow(sub_matrix))


# Apply t-test and linear model for each pathway (row)
for (i in 1:nrow(sub_matrix)) {
  # Get z-scores for this pathway across samples
  a_values <- sub_matrix[i, ,drop = FALSE]

  mean_total[i] <- mean(a_values, na.rm = TRUE)
  
  # Compute standard deviation for each group
  std_total[i] <- sd(a_values, na.rm = TRUE)
 
  
  if (i %% 10 == 0) {
    cat(blue(i), "out of", green(nrow(sub_matrix)), " ", red(round(i/nrow(sub_matrix)*100, 2)), "%", "\n")
  }
  
  
}



# Create results tibble
results <- tibble(
  Pathway = rownames(sub_matrix),
  mean_total=mean_total,
  std_total = std_total
)




# Create the output list with metadata
analysis_output <- list(
  results = results,
  metadata = list(
    date = Sys.time(),
    group="all samples",
    group_samples = colnames(sub_matrix),

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



saveRDS(analysis_output, paste0("nf1g/gsvas/ds/hu-gsva_pathway_stats_",
        
        "all_samples-minus_wt", 
        
        ".rds")

)



