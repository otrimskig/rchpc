source("libs.R")
library(tidyverse)
library(dtplyr)


all_counts_tidy<-readRDS("k19mf/ds/vm-01-gene_id_hu_m_rpkms.rds")


sample_info<-readRDS("k19mf/ds/vm-00-sample_info.rds")%>%
  select(sample_id, mouse_num, tumor_type)


file<-all_counts_tidy%>%
  left_join(sample_info)%>%
  mutate(column_name=paste0(tumor_type, "_", mouse_num, "_", sample_id))%>%
  
  select(-sample_id, -mouse_num, -tumor_type)%>%
  
  filter(!is.na(gene_name_ms))%>%
  
  pivot_wider(names_from=column_name, values_from = c(read_count, rpkm),  names_glue = "{column_name}_{.value}")%>%
  
  select(1:4, sort(names(.)[5:last_col()]))



write_csv(file, "k19mf/ds/foot_pads_norm_counts.csv", na = "")
