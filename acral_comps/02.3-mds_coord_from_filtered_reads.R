source("libs.R")
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(edgeR)


vm02<-readRDS("acral_comps/ds/vm-02-filtered_rpkms.rds")

sample_info<-readRDS("acral_comps/ds/v00-sample_info.rds")%>%
  mutate(exclude=if_else(sample_id=="x"))


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
saveRDS(p, "acral_comps/ds/mds_dge.rds")



