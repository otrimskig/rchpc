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


heatmap(as.matrix(hm))






de <- read_csv("de-nf1 wt; pten wt; ink wt; atrx wt vs. nf1 KO; pten KO; ink KO; atrx KO.csv")

hm<-de%>%
  filter(!is.na(gene_name))%>%
  
  filter(abs(logFC)>=2)%>%
  
  
  filter(FDR<.001)%>%
  
  
  sample_n(500)%>%
  
  
  
  group_by(gene_name)%>%
  slice(1)%>%
  ungroup()%>%
  
  
  
  
  column_to_rownames("gene_name")%>%
  
  select(starts_with("rpkm"))




  





sample_info<-readRDS("v05-all_counts_plus_info.rds")%>%
  group_by(sample_id)%>%
  slice(1)%>%
  ungroup()%>%
  select(sample_id, mouse_num:patho_cat)

sample_info$resultant_geno

ha<-HeatmapAnnotation(df=sample_info)





ha <- HeatmapAnnotation(foo = 1:13, col = list(foo = col_fun))






h<-Heatmap(as.matrix(hm), show_row_dend = FALSE, row_dend_reorder = T)






ha2<- HeatmapAnnotation(bar = sample(letters[1:3], 10, replace = TRUE))


hm<-de%>%
  filter(!is.na(gene_name))%>%
  
  filter(abs(logFC)>=1)%>%
  
  
  filter(FDR<.001)%>%
  
  sample_n(100)%>%
  
  
  group_by(gene_name)%>%
  slice(1)%>%
  ungroup()%>%
  
  
  
  
  column_to_rownames("gene_name")%>%
  
  select(starts_with("rpkm"))


heatmap(as.matrix(hm))



library(ggplot2)
library(hrbrthemes)

hm%>%
  rownames_to_column(var="gene")%>%
  
  #sample_n(100)%>%
  
  
  pivot_longer(2:last_col(), names_to = "sample_id", values_to = "rpkm")%>%
  
  
  mutate(sample_id=sub("^rpkm_", "", sample_id))%>%
  
  
  
  ggplot(aes(x=sample_id, y=gene, fill=rpkm))+
  geom_tile()+
  scale_fill_viridis(discrete=FALSE)



head()
