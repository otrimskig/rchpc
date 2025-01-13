source("libs.R")
library(tidyverse)
library(edgeR)
library(dtplyr)

sample_info<-readRDS("acral_paired/ds/v00-sample_info.rds")
read_counts<-readRDS("acral_paired/ds/v02-filtered_rpkms.rds")


name_map<-readRDS("acral_paired/ds/name_map.rds")%>%
  select(sample_id, mouse_num_type_clean)


all_data<-read_counts%>%left_join(sample_info)%>%left_join(name_map)%>%
  ungroup()%>%arrange(gene_name_ms)%>%
  select(gene_id_ms, gene_name_ms, gene_name_hu, read_count, rpkm, mouse_num_type_clean)%>%
  
  pivot_wider(values_from = c(read_count, rpkm), names_from = mouse_num_type_clean)%>%
  arrange(gene_name_ms)

openxlsx::write.xlsx(all_data, "acral_paired/outs/acral_paired_rawc_rpkm_filtered1_wide.xlsx")
