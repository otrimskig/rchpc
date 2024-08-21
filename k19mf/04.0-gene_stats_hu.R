source("libs.R")

library(tidyverse)
library(dtplyr)


data<-readRDS("k19mf/ds/vm-02-filtered_rpkms.rds")




#make dataset file for unique human genes, and their stats.
#keep highest expressed version (by mean exp). 
#keep link to ms genes in form of gene ids and names.
#you can use this if nec to link back to previous expression data with ms. 
#(not that that should be nec.)

unique_hu_genes<-readRDS("k19mf/ds/gene_stats.rds")%>%
  filter(!is.na(gene_name_hu))%>%
  group_by(gene_name_hu)%>%
  arrange(desc(mean_rpkm))%>%
  slice(1)%>%ungroup()


saveRDS(unique_hu_genes, "k19mf/ds/unique_hu_gene_stats.rds")




#make a wide-version df of rpkms for unique human genes, as 
#filtered by the above unique list of genes.
#column titles are mouse numbers. 

mouse_nums<-readRDS("k19mf/ds/vm-00-sample_info.rds")%>%
  select(mouse_num, sample_id)

data_hu_only<-data%>%
  semi_join(unique_hu_genes, by="gene_id_ms")%>%
  
  left_join(mouse_nums)%>%
  select(-sample_id)%>%
  select(gene_name_hu, rpkm, mouse_num)%>%
  filter(!is.na(gene_name_hu))%>%
  
  pivot_wider(names_from = mouse_num, values_from=rpkm)



saveRDS(data_hu_only, "k19mf/ds/vm-h-01-rpkms_wide_human.rds")





