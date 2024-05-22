library(tidyverse)
library(ggplot2)
library(ggrepel)
#library(viridis)
library(edgeR)


vm02<-readRDS("ds/vm-01-gene_id_hu_m_rpkms.rds")

vm02<-vm02%>%
  filter(sample_id!="X23908X21")%>%
  filter(sample_id!="X23908X15")%>%
  filter(sample_id!="X23908X07")



sample_info<-readRDS("ds/v07-per_sample_info.rds")


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




