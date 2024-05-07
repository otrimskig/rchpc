library(tidyverse)


setwd("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig")




#read in all data. tidy dataframe, with all sample info
#plus every raw count for each sample and gene.
#remove 15 and 21 for low quality. 
all_data<-readRDS("23908R/v05-all_counts_plus_info.rds")%>%
  
  #rm x15 and x21 from dataset, the 2 samples that performed poorly in QC. 
  filter(!str_ends(sample_id, "x15|x21"))%>%
  rename(gene_name=external_gene_name)%>%
  
  filter(!is.na(gene_name))%>%
  
  group_by(sample_id, gene_name)%>%slice(1)%>%ungroup()%>%
  group_by(sample_id, gene_id)%>%slice(1)%>%ungroup()%>%
  
  mutate(patho_cat2=patho_cat)%>%
  mutate(patho_cat=if_else(grepl("^[0-9]", patho_cat), substr(patho_cat, 1,1), patho_cat))%>%
  
  mutate(patho_grade=if_else(grepl("^[0-9]", patho_cat2), substr(patho_cat2, 3,3), patho_cat2))


saveRDS(all_data, "23908R/v06-all_counts_plus_info.rds")
