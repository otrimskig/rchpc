library(tidyverse)
library(dtplyr)

v6_ms<-readRDS("ds/v06-all_counts_plus_info.rds")%>%
  select(gene_name)%>%
  unique()%>%
  filter(!is.na(gene_name))%>%
  rename(gene_name_ms=gene_name)



library(biomaRt)



human<- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org")

mouse<-useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org")

gene_conv<- getLDS(attributes = c("hgnc_symbol"), 
                 filters = "hgnc_symbol", 
                 values = v6_ms$gene_name_ms , 
                 mart = human, 
                 attributesL = c("mgi_symbol"), 
                 martL = mouse, uniqueRows=T)



v6_mh<-v6_ms%>%
  full_join(gene_conv%>%rename(gene_name_ms=MGI.symbol, gene_name_hu=HGNC.symbol))%>%
  filter(!is.na(gene_name_hu))
  



v8<-readRDS("ds/v06-all_counts_plus_info.rds")%>%
  rename(gene_name_ms=gene_name)%>%
  left_join(v6_mh)%>%
  relocate(gene_name_hu, .after=gene_name_ms)




saveRDS(v8, "ds/v08-rpkms_w_meta_ms_hu.rds")





