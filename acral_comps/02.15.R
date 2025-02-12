source("libs.R")
library(tidyverse)
library(dtplyr)


#relative to rchpc proj
proj_dir<-"acral_comps/"


df1<-readRDS("acral_comps/ds/v01-gene_id_rpkms.rds")

df2<-readRDS("acral_paired/ds/v01-gene_id_hu_m_rpkms.rds")%>%
  ungroup()%>%
  dplyr::select(c(1:4))%>%
  unique()


df3<-df2%>%left_join(df1)

saveRDS(df3, "acral_comps/ds/vm-01-gene_id_hu_m_rpkms.rds")
