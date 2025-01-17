##########################
source("libs.R")
library(tidyverse)
library(edgeR)
library(dtplyr)



df0<-readRDS("acral_paired/ds/v01-gene_id_hu_m_rpkms.rds")


df1<-df0%>%
  pivot_wider(names_from = sample_id, values_from = c(read_count, rpkm))




write_csv(df1, "acral_paired/outs/counts_matrix_raw_plus_rpkms.csv", na = "")



