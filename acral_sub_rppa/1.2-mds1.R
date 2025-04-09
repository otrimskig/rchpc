source("libs.R")
library(ggplot2)
library(ggrepel)
#library(viridis)
library(edgeR)




vm0<-readRDS("acral_sub_rppa/ds/rppa_df_list0.rds")

sample_info0<-readRDS("acral_sub_rppa/ds/sample_info0.rds")


sample_interest<-sample_info0%>%
  filter(sample_type!="none")

vm1<-vm0[["L4 (linear)"]][["df_long"]]


mat<-vm1%>%
  filter(order %in% sample_interest$order)%>%
  
  filter(order!=455&order!=456)%>%
  filter(order!=469)%>%
  
  pivot_wider(names_from = antibody_name, values_from=rppa_value)%>%
  column_to_rownames("order")%>%
  as.matrix.data.frame()%>%
  t()

mat_num <- apply(mat, 2, as.numeric)




p<-plotMDS(mat_num, labels=colnames(mat))


saveRDS(p, "acral_sub_rppa/ds/mds_dge0.rds")



