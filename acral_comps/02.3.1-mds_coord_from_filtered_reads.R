source("libs.R")
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(edgeR)


vm01<-readRDS("acral_comps/ds/vm-02-filtered_rpkms.rds")

sample_info<-readRDS("acral_comps/ds/v00-sample_info.rds")%>%
  mutate(exclude=NA)%>%
  
  # mutate(exclude=if_else(sample_id=="x14154x08"|
  #                          sample_id=="x14154x07", 1, NA))%>%
  
  mutate(exclude=if_else(assigned_percent<20, 1, NA))%>%
  filter(is.na(exclude))

vm02<-vm01%>%semi_join(sample_info, by="sample_id")



mat<-vm02%>%
  select(gene_id_ms, sample_id, rpkm)%>%
  pivot_wider(names_from = sample_id, values_from = rpkm)%>%
  column_to_rownames("gene_id_ms")%>%
  as.matrix.data.frame()

lib_size<-vm02%>%
  select(sample_id, read_count)%>%
  ungroup()%>%
  group_by(sample_id)%>%
  summarise(lib_size=sum(read_count))%>%
  ungroup()%>%
  arrange(sample_id)%>%
  pull(lib_size)

dge<-DGEList(counts = mat, 
        lib.size = lib_size)



p<-plotMDS.DGEList(dge)

plotMDS.DGEList(dge)

stop()
saveRDS(p, "acral_comps/ds/mds_dge05.rds")



