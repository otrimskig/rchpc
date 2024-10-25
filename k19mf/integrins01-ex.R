library(tidyverse)
library(dtplyr)
library(broom)

df<-readRDS("k19mf/ds/vm-02-filtered_rpkms.rds")

gs<-readRDS("k19mf/ds/gene_stats.rds")


its<-gs%>%
  filter(grepl("^Itg", gene_name_ms))%>%
  pull(gene_name_ms)




its_df<-df%>%filter(gene_name_ms %in% its)
  


# its_df%>%
#   group_by(gene_name_ms)%>%count()%>%view()




saveRDS(its_df,"k19mf/ds/intds.rds")



tumor_type<-readRDS("k19mf/ds/vm-00-sample_info.rds")%>%select(sample_id, mouse_num, tumor_type)

its_ps<-its_df%>%select(gene_name_ms, sample_id, rpkm)%>%
  left_join(tumor_type)%>%
  mutate(tumor_type=gsub("foot pad", "acral", tumor_type))%>%
  
  filter(mouse_num!="3126")


op1<-its_ps%>%
  filter(gene_name_ms=="Itga1")


op1

results <- its_ps %>%
  group_by(gene_name_ms) %>%

  summarise(t_test = list(t.test(rpkm ~ tumor_type))) %>%
  mutate(p_value = map_dbl(t_test, ~ .x$p.value)) %>%
  mutate(p_adj_bonferroni = p.adjust(p_value, method = "bonferroni"),  # Bonferroni correction
         p_adj_fdr = p.adjust(p_value, method = "fdr"))%>%
  select(-t_test)
  #select(gene_name_ms, p_value)


saveRDS(its_ps%>%
  left_join(results),
  
  "k19mf/ds/int-pvals-ex.rds")
