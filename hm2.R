library(tidyverse)

library(ComplexHeatmap)











de <- read_csv("de-nf1 KO; pten KO; ink KO; atrx wt vs. nf1 KO; pten KO; ink KO; atrx KO.csv")

hm<-de%>%
  filter(!is.na(gene_name))%>%
  
  filter(abs(logFC)>=1.5)%>%
  
  
  filter(FDR<.1)%>%
  
  group_by(gene_name)%>%
  slice(1)%>%
  ungroup()%>%  
  column_to_rownames("gene_name")%>%
  
  select(starts_with("rpkm"))




hp<-Heatmap(as.matrix(hm))
  
sample_info2<-sample_info%>%
  mutate(col_name=paste0("rpkm_", sample_id))





all<-readRDS("v05-all_counts_plus_info.rds")%>%
  
  rename(gene_name=external_gene_name)%>%
  filter(!is.na(gene_name))%>%
  
  
  #filter(abs(logFC)>=1.5)%>%
  
  
  #filter(FDR<.1)%>%
  
  group_by(mouse_num, gene_name)%>%
  slice(1)%>%
  ungroup()
  

all2<-all%>%
  select(gene_name, mouse_num, rpkm)%>%
  
  pivot_wider(values_from = rpkm, names_from = mouse_num)
    
    
    
all_hm<-Heatmap(as.matrix(all2))




