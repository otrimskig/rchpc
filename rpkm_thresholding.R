library(tidyverse)
library(dtplyr)

#set filter to remove bad samples.
remove<-as.character(c("x23908x07","x23908x15","x23908x21"))


#get all reads for all samples.
df<-readRDS("ds/vm-01-gene_id_hu_m_rpkms.rds")%>%
  mutate(sample_id=tolower(sample_id))%>%
  filter(!sample_id %in% remove)%>%
  
  ungroup()%>%
  
  #calculate some stats for each gene to filter.
  group_by(gene_id_ms)%>%
  mutate(mean_rpkm = mean(rpkm))%>%
  mutate(max_rpkm = max(rpkm))%>%
  ungroup()%>%
  
  #remove genes where the max rpkm for any sample was at least 1.
  filter(max_rpkm>=1)
  
  
  
#to be used for de expression sets.
#remove any unnamed ms genes.
#ensure only unique mouse gene names are kept.
#but keep the one with the highest average rpkm.
df2<-df%>%
  filter(!is.na(gene_name_ms))%>%
  group_by(sample_id, gene_name_ms)%>%
  arrange(desc(mean_rpkm))%>%
  slice(1)%>%
  select(-mean_rpkm, -max_rpkm)%>%
  ungroup()

saveRDS(df2, "ds/vm-02-filtered_rpkms.rds")




#get stats for all gene name - id combos. Some have no mouse name, but have a human name. 
#some are duplicated. Keep a dataset with all unique combos with rpkm stats for each.
gene_stats<-df%>%
  select(gene_id_ms, gene_name_hu, gene_name_ms, mean_rpkm, max_rpkm)%>%
  unique()

saveRDS(gene_stats, "ds/gene_stats.rds")
