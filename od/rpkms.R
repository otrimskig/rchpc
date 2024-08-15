library(tidyverse)
        

setwd("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/23908R")

df<-readRDS("v02-counts_with_geneIDS.rds")






df2<-df%>%
  pivot_longer(cols = 5:last_col(), 
               names_to="sample_name", 
               values_to = "raw_read_count")%>%
  rename(gene_length=length)


total_reads<-df2%>%
  group_by(sample_name)%>%
  summarise(total_reads=sum(raw_read_count))


df3<-df2%>%
  left_join(total_reads)%>%
  
  mutate(rpkm = raw_read_count/(gene_length/1000*total_reads/1E6))



genes_above_1<-df3%>%
  group_by(gene_id)%>%
  summarise(max_rpkm=max(rpkm))%>%
  filter(max_rpkm>1)



df4<-df3%>%
  semi_join(genes_above_1)



saveRDS(df4, "v03-rpkms.rds")




max_dif<-df4%>%
  group_by(gene_id, external_gene_name)%>%
  summarise(rpkm_dif=max(rpkm)-min(rpkm))%>%
  arrange(desc(rpkm_dif))%>%
  filter(!is.na(external_gene_name))%>%
  ungroup()%>%
  mutate(name = fct_reorder(external_gene_name, desc(rpkm_dif)))





