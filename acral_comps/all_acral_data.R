source("libs.R")
library(tidyverse)


df_k19mf<-readRDS("k19mf/ds/vm-02-filtered_rpkms.rds")
sa_k19mf<-readRDS("k19mf/ds/vm-00-sample_info.rds")%>%
  rename(sample_type=tumor_type)%>%
  relocate(sample_id, mouse_num, sample_type, genotype, contains("reads"), contains("assign"))


df_acpair<-readRDS("acral_paired/ds/v02-filtered_rpkms.rds")
sa_acpair<-readRDS("acral_paired/ds/v00-sample_info.rds")%>%
  mutate(genotype=NA)%>%
  relocate(sample_id, mouse_num, sample_type, genotype, contains("reads"), contains("assign"))




sa_integrated<-full_join(sa_acpair, sa_k19mf)

saveRDS(sa_integrated, "acral_comps/ds/v00-sample_info.rds")




g01<-df_k19mf%>%select(gene_id_ms)%>%unique()
g02<-df_acpair%>%select(gene_id_ms)%>%unique()




counts1<-readRDS("k19mf/ds/vm-00-all_counts_sample_names_fixed.rds")
counts2<-readRDS("acral_paired/ds/v00-all_counts.rds")%>%
  rename(gene_id=GeneID, length=Length)

counts_int<-left_join(counts1, counts2)

saveRDS(counts_int, "acral_comps/ds/v00-all_counts.rds")
