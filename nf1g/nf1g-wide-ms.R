source("libs.R")

library(tidyverse)
library(dtplyr)


df<-readRDS("nf1g/ds/vm-01-gene_id_hu_m_rpkms.rds")


df2<-df%>%
  select(-c("gene_len_ms", "gene_name_ms", "gene_name_hu", "read_count"))%>%
  pivot_wider(names_from="sample_id", values_from = "rpkm")%>%
  janitor::clean_names()%>%
  
  column_to_rownames("gene_id_ms")%>%
  as.matrix.data.frame()


saveRDS(df2, "nf1g/ds/ms_ens_rpkms_wide.rds")
