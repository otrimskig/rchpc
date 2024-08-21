source("libs.R")

#cudc rna seq timex analysis 
#timex does not care about comparisons, so any group loaded is fine 
#needs RPKM data in a format where the rownames are the gene names (gene symbol) and the rpkms are a matrix 


library(tidyverse)
library(dtplyr)
library(GSVA)



data<-readRDS("ds/vm-02-filtered_rpkms.rds")

unique_hu_genes<-readRDS("ds/gene_stats.rds")%>%
  filter(!is.na(gene_name_hu))%>%
  group_by(gene_name_hu)%>%
  arrange(desc(mean_rpkm))%>%
  slice(1)%>%ungroup()


mouse_nums<-readRDS("ds/v10-per_sample_updated.rds")%>%
  select(mouse_num, sample_id)


data_hu_only<-data%>%
  semi_join(gene_stats, by="gene_id_ms")%>%
  left_join(mouse_nums)%>%
  select(-sample_id)%>%
  select(gene_name_hu, rpkm, mouse_num)%>%
  pivot_wider(names_from = mouse_num, values_from=rpkm)


saveRDS(data_hu_only, "ds/vm-h-01-rpkms_wide_human.rds")

mat<-data_hu_only%>%
  column_to_rownames("gene_name_hu")%>%
  data.matrix()










load("timex/ds/allSignatures.rda")  


onco<-qusage::read.gmt("timex/ds/c6.all.v2024.1.Hs.symbols.gmt")

names(onco) <- paste0("onco_", names(onco))


Signature_list <- c(Hallmark, kegg, TIMEx, Immune_sig, onco)

