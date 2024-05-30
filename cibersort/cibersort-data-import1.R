library(biomaRt)
library(dtplyr)
library(tidyverse)

#all genes that appear in reference dataset. 
ref_genes<-data.table::fread("cibersort/41598_2017_BFsrep40508_MOESM313_ESM.txt")%>%
  as_tibble()%>%
  select(1)%>%
  rename(gene_name_ms=V1)%>%
  mutate(ref_gene_name_ms=gene_name_ms)%>%
  mutate(ref_gene_id=1:n())


#list of genes that we have in RNA dataset.
#readRDS("ds/vm-01-gene_id_hu_m_rpkms.rds")

seq_genes<-readRDS("ds/vm-01-gene_id_hu_m_rpkms.rds")%>%filter(sample_id=="X23908X01")%>%
  select(-sample_id, -read_count, -rpkm)%>%
  group_by(gene_id_ms)%>%slice(1)%>%ungroup()%>%
  mutate(seq_gene_name_ms=gene_name_ms)%>%
  mutate(seq_exp_id=1:n())





mouse<-useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#biomart results
ref_genes_bm_r<-getBM(attributes = c("entrezgene_id", "external_gene_name",
                                     "external_synonym", "ensembl_gene_id"),
                     
                     values = ref_genes$gene_name_ms, 
                     mart = mouse)%>%
  
  rename(ref_gene_name_ms=external_gene_name)%>%
  mutate(across(everything(), ~ na_if(as.character(.), "")))



ref_genes_bm<-ref_genes%>%left_join(ref_genes_bm_r)%>%
  pivot_longer(cols=c("gene_name_ms", "external_synonym"), values_to="gene_name_ms", names_to = "type")
  


#biomart results
seq_genes_bm_r<-getBM(attributes = c("entrezgene_id", "external_gene_name",
                                   "external_synonym", "ensembl_gene_id"),
                    
                    values = seq_genes$gene_id_ms, 
                    mart = mouse)%>%
  
  rename(seq_gene_name_ms=external_gene_name)%>%
  mutate(across(everything(), ~ na_if(as.character(.), "")))


seq_genes_bm<-seq_genes%>%left_join(seq_genes_bm_r, relationship = "many-to-many")%>%
  pivot_longer(cols=c("gene_name_ms", "external_synonym"), values_to="gene_name_ms", names_to = "type")
  




match_check<-ref_genes_bm%>%
  mutate(ref_gene="1")%>%
  full_join(seq_genes_bm%>%select(type, gene_name_ms, seq_exp_id))%>%
  filter(ref_gene=="1")%>%
  full_join(seq_genes_bm%>%select(entrezgene_id, seq_exp_id))%>%
  filter(ref_gene=="1")%>%
  full_join(seq_genes_bm%>%select(ensembl_gene_id, seq_exp_id))%>%
  filter(ref_gene=="1")%>%
  
  select(ref_gene_id, seq_exp_id)%>%
  unique()








match1<-ref_genes_bm%>%left_join(seq_genes, by="gene_name_ms")%>%
  select(-gene_name_ms)


















all_genes_seq<-unfiltered_reads%>%filter(sample_id=="X23908X01")%>%select(-sample_id, -read_count, -rpkm)%>%
  mutate(seq_exp_id=gene_id_ms)%>%ungroup()




ref_genes%>%semi_join(all_genes_seq, by="gene_name_ms")%>%
  group_by(gene_name_ms)%>%slice(1)%>%view()

ref_genes%>%anti_join(all_genes_seq, by="gene_name_ms")%>%
  group_by(gene_name_ms)%>%slice(1)%>%view()







df1<-ref_genes%>%semi_join(all_genes_seq, by="gene_name_ms")%>%
  left_join(all_genes_seq, by="gene_name_ms")



no_gene_name_match<-ref_genes%>%anti_join(all_genes_seq, by="gene_name_ms")














all_genes<-left_join(all_genes_seq, all_genes_bm)

# listAttributes(mouse)%>%view()

ref_genes_biomart<-left_join(ref_genes, ref_gene_info)

lost_genes<-ref_genes_biomart%>%
  filter(is.na(entrezgene_id))%>%
  filter(is.na(ensembl_gene_id))%>%
  mutate(fake_human=toupper(gene_name_ms))









human<-useMart("ensembl", dataset = "hsapiens_gene_ensembl")



poop<-getBM(attributes = c("entrezgene_id", "external_gene_name",	
                     "ensembl_transcript_id", "external_synonym", "ensembl_gene_id", "hgnc_id"),
      values = lost_genes$fake_human, 
      mart = human)
















df2<-no_gene_match_alts3%>%semi_join(all_genes_seq, by="gene_name_ms")%>%
  left_join(all_genes_seq, by="gene_name_ms")%>%
  rename(ref_ms=external_synonym)%>%
  select(ref_ms, seq_exp_id)




all_genes_seq%>%semi_join(df1,by="seq_exp_id")%>%group_by(seq_exp_id)%>%slice(1)%>%ungroup()->match1

all_genes_seq%>%semi_join(df2,by="seq_exp_id")%>%group_by(seq_exp_id)%>%slice(1)%>%ungroup()->match2


full_join(match1, match2)%>%group_by(seq_exp_id)%>%slice(1)%>%ungroup()->all_matched






unmatched58<-ref_genes%>%
  anti_join(matched)



no_gene_match_alts2<-no_gene_match_alts1%>%
  count(external_gene_name)



df2<-df1%>%filter(is.na(seq_exp_id))%>%select(-seq_exp_id)%>%
  rename(ref_ms=gene_name_ms)%>%
  rename(gene_name_ms=gene_name_ms_syn)

df3<-df2%>%
  filter(!is.na(gene_name_ms))%>%
  left_join(all_genes_seq%>%select(gene_name_ms, seq_exp_id))%>%
  filter(!is.na(seq_exp_id))




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
