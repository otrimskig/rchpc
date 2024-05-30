library(biomaRt)
library(dtplyr)
library(tidyverse)

#all genes that appear in reference dataset. 
ref_genes<-data.table::fread("cibersort/41598_2017_BFsrep40508_MOESM313_ESM.txt")%>%
  as_tibble()%>%
  select(1)%>%
  rename(gene_name_ms=V1)

#list of genes that we have in RNA dataset.
unfiltered_reads<-readRDS("ds/vm-01-gene_id_hu_m_rpkms.rds")

all_genes_seq<-unfiltered_reads%>%filter(sample_id=="X23908X01")%>%select(-sample_id, -read_count, -rpkm)%>%
  mutate(seq_exp_id=gene_id_ms)%>%ungroup()


mouse<-useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org")

ref_gene_info<-getBM(attributes = c("entrezgene_id", "external_gene_name", "external_synonym"),
                 filters="external_gene_name",
                 values = ref_genes$gene_name_ms, 
                 mart = mouse)%>%
  rename(gene_id_ms=entrezgene_id,
         gene_name_ms=external_gene_name,
         gene_name_ms_syn=external_synonym)%>%
  mutate(gene_name_ms_syn=if_else(gene_name_ms_syn=="", NA, gene_name_ms_syn))




#get any genes from ref genes that didn't have hits in biomart.
no_search_result<-ref_gene_info%>%
  anti_join(ref_genes, by="gene_name_ms")

#should be TRUE
count(no_search_result)==0





df1<-ref_gene_info%>%
  left_join(all_genes_seq)

df2<-df1%>%filter(is.na(seq_exp_id))%>%select(-seq_exp_id)%>%
  rename(ref_ms=gene_name_ms)%>%
  rename(gene_name_ms=gene_name_ms_syn)

df3<-df2%>%
  left_join(all_genes_seq%>%select(gene_name_ms, seq_exp_id))





# 
# all_genes_seq%>%
#   semi_join(ref_gene_info, by="gene_name_ms")
# 
# unmatched_genes<-ref_gene_info%>%
#   anti_join(all_genes_seq, by="gene_name_ms")%>%
#   group_by(gene_name_ms)%>%slice(1)
# 
# 
# 
# 
# unmatched_genes2<-ref_gene_info%>%
#   anti_join(all_genes_seq, by="gene_id_ms")%>%
#   group_by(gene_name_ms)%>%slice(1)
# 
# 
# unmatched_genes3<-ref_gene_info%>%
#   semi_join(all_genes_seq, by=c("gene_name_ms_syn"="gene_name_ms"))%>%
#   group_by(gene_name_ms)%>%slice(1)
# 
# 




stop("done bit")


# 
# human<- useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org")
# 
# 
# hu_gene_info<- getLDS(attributes = c("hgnc_symbol"), 
#                    filters = "hgnc_symbol", 
#                    values = ref_genes$gene_name_ms, 
#                    mart = human, 
#                    attributesL = c("mgi_symbol"), 
#                    martL = mouse, uniqueRows=T)
