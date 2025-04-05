source("libs.R")
# suppressPackageStartupMessages(library(ComplexHeatmap))
# ht_opt$message = FALSE
library(gridtext)
library(circlize)
library(purrr)
library(tidyverse)
library(ggplot2)



sample_info<-readRDS("nf1g/ds/v10-per_sample_updated.rds")
matu0<-readRDS("timex/ds/hu-msig_all-gsva-values.rds")

gsva_comp_file_path<-"nf1g/gsvas/ds/hu-gsva_pathway_stats_a86_Gliomas][nf1 KO; pten KO; ink KO; atrx KOvsnf1 KO; pten KO; ink KO; atrx wt.rds"



gsva_comp<-readRDS(gsva_comp_file_path)

sa_v<-c(gsva_comp$metadata$group_a_samples,
  gsva_comp$metadata$group_b_samples
)

sa_mat0<-matu0[,sa_v, drop=FALSE]

if(exists("results2", where = gsva_comp)){
  
  comp_mat<-gsva_comp$results2
  
}else{
  
  comp_mat<-gsva_comp$results
} #define comp_mat based on where results were stored in list.


group_a_name<-gsva_comp$metadata$group_a
group_b_name<-gsva_comp$metadata$group_b

subset_name <- sub("\\].*", "", 
              sub("a86_", "", 
                  sub("hu-gsva_pathway_stats_", "", basename(gsva_comp_file_path))
              ))


# Assuming your matrix is named `my_matrix` and it has columns 'mean_group_a' and 'mean_group_b'
n <- 50  # Specify the number of top rows you want to extract

# Get top n rows for each group, preserving row indices
top_a_indices <- order(comp_mat$mean_group_a, decreasing = TRUE)[1:n]
top_b_indices <- order(comp_mat$mean_group_b, decreasing = TRUE)[1:n]

# Combine indices and remove duplicates (if any)
top_indices <- unique(c(top_a_indices, top_b_indices))

# Subset the matrix using these indices
top_hits <- comp_mat[top_indices, ]%>%
  as_tibble()




sa_df1 <- sa_mat0[top_hits$Pathway, , drop = FALSE] %>%
  as_tibble(rownames = "Pathway")%>%
  pivot_longer(cols= -Pathway, names_to = "mouse_num", values_to = "z_score")%>%
  left_join(sample_info%>%
              select(mouse_num, patho_cat_name, resultant_geno))%>%
  left_join(top_hits%>%select(Pathway, mean_total))


sa_df2<-sa_df1

ggplot(data=sa_df2, aes(x=z_score, y=reorder(Pathway, mean_total)))+
  geom_point(aes(color=resultant_geno), size=3, alpha=.5)+
  theme_bw()+
  facet_grid(cols = vars(resultant_geno)) 



##################



# Subset the matrix using these indices
top_diff <- comp_mat%>%
  as_tibble()%>%
  #filter(min_pval<.05)%>%
  arrange(desc(abs_diff))%>%
  slice(1:50)




diff_df1 <- sa_mat0[top_diff$Pathway, , drop = FALSE] %>%
  as_tibble(rownames = "Pathway")%>%
  pivot_longer(cols= -Pathway, names_to = "mouse_num", values_to = "z_score")%>%
  left_join(sample_info%>%
              select(mouse_num, patho_cat_name, resultant_geno))
  


diff_df2<-diff_df1%>%
  left_join(top_diff%>%select(Pathway, abs_diff, mean_total, diff))

ggplot(data=diff_df2, aes(x=z_score, y=reorder(Pathway, diff)))+
  geom_point(aes(color=resultant_geno), size=3, alpha=.5)+
  theme_bw()
