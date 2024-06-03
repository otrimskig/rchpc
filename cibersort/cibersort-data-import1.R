library(biomaRt)
library(tidyverse)
library(dtplyr)
library(fuzzyjoin)

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
  pivot_longer(cols=c("gene_name_ms", "external_synonym"), values_to="gene_name_ms", names_to = "type1")%>%
  mutate(gene_name_ms=tolower(gene_name_ms))
  


#biomart results
seq_genes_bm_r<-getBM(attributes = c("entrezgene_id", "external_gene_name",
                                   "external_synonym", "ensembl_gene_id"),
                    
                    values = seq_genes$gene_id_ms, 
                    mart = mouse)%>%
  
  rename(seq_gene_name_ms=external_gene_name)%>%
  mutate(across(everything(), ~ na_if(as.character(.), "")))


seq_genes_bm<-seq_genes%>%left_join(seq_genes_bm_r, relationship = "many-to-many")%>%
  pivot_longer(cols=c("gene_name_ms", "external_synonym"), values_to="gene_name_ms", names_to = "type2")%>%
  mutate(gene_name_ms=tolower(gene_name_ms))
  


seq_gene_name_ms<-seq_genes_bm%>%select(gene_name_ms, seq_exp_id)%>%unique()%>%
  filter(!is.na(gene_name_ms))%>%semi_join(ref_genes_bm)

seq_entrez<-seq_genes_bm%>%select(entrezgene_id, seq_exp_id)%>%unique()%>%
  filter(!is.na(entrezgene_id))%>%semi_join(ref_genes_bm)

seq_ensembl<-seq_genes_bm%>%select(ensembl_gene_id, seq_exp_id)%>%
  filter(!is.na(ensembl_gene_id))%>%unique()%>%semi_join(ref_genes_bm)




match_check<-ref_genes_bm%>%
  mutate(ref_gene="1")%>%
  left_join(seq_gene_name_ms)%>%
  left_join(seq_entrez)%>%
  left_join(seq_ensembl)%>%
  filter(ref_gene=="1")%>%
  

  
  #select(ref_gene_id, seq_exp_id)%>%
  unique()



matches<-match_check%>%
  group_by(ref_gene_id)%>%
  mutate(match=if_else(any(!is.na(seq_exp_id)), "1", "0"))



matches2<-matches%>%filter(match=="1")

unmatched<-matches%>%filter(match=="0")%>%
  filter(!is.na(gene_name_ms))%>%
  select(ref_gene_id, ref_gene_name_ms)%>%
  unique()

library(googlesheets4)
gs4_auth(email = "gotrimski@gmail.com")
unmatched%>%
    range_write("https://docs.google.com/spreadsheets/d/1A4wa8WsbsazEfBSndcSqgEk1BMZokdTpIPNKjHV5PWE/edit#gid=378460454",
                sheet = "gns",
                .,
                 reformat=FALSE,
                range = "A1")




#manual annotation. 
gs4_auth(email = "gotrimski@gmail.com")
unmatched_manual_edits<-read_sheet("https://docs.google.com/spreadsheets/d/1A4wa8WsbsazEfBSndcSqgEk1BMZokdTpIPNKjHV5PWE/edit#gid=378460454",
           sheet="gns")%>%
 
  mutate(across(everything(), ~ na_if(as.character(.), "")))%>%
  mutate(across(1:last_col(), function(x){na_if(x, "NULL")}))

unmatched_ids<-unmatched_manual_edits%>%
  mutate(mgi_id=if_else(!is.na(MGI),paste0("MGI:", MGI), NA))%>%

  select(-MGI)%>%
  rename(external_synonym=alias,
         ensembl_gene_id=ens)



#biomart results
seq_genes_bm_r2<-getBM(attributes = c("entrezgene_id", "mgi_id","external_gene_name",
                                     "external_synonym", "ensembl_gene_id", "hgnc_id"),
                      
                      values = seq_genes$gene_id_ms, 
                      mart = mouse)%>%
  
  #rename(seq_gene_name_ms=external_gene_name)%>%
  
  mutate(across(everything(), ~ na_if(as.character(.), "")))%>%
  
  
  full_join(seq_genes%>%mutate(entrezgene_id=as.character(gene_id_ms)), by="entrezgene_id")



#now check against all unmatched results for like ids. 
mgi_matches<-unmatched_ids%>%
  left_join(seq_genes_bm_r2, na_matches="never", by="mgi_id")%>%
  filter(!is.na(seq_exp_id))


syn_matches_hu<-unmatched_ids%>%
  filter(is.na(mgi_id))%>%
  left_join(seq_genes_bm_r2, na_matches ="never", by=c("external_synonym"="gene_name_hu"))%>%
  filter(!is.na(seq_exp_id))



name_matches<-unmatched_ids%>%
  filter(is.na(mgi_id))%>%
  filter(is.na(external_synonym))%>%
  left_join(seq_genes_bm_r2, na_matches ="never", by=c("ref_gene_name_ms"="seq_gene_name_ms"))%>%
  filter(!is.na(seq_exp_id))






unmatched_bm2<-getBM(attributes = c("entrezgene_id", "mgi_id","external_gene_name",
                                      "external_synonym", "ensembl_gene_id", "hgnc_id"),
                       
                       values = unmatched_ids%>%select(-ref_gene_id), 
                       mart = mouse)%>%
  
  #rename(seq_gene_name_ms=external_gene_name)%>%
  
  mutate(across(everything(), ~ na_if(as.character(.), "")))
  
  
  full_join(seq_genes%>%mutate(entrezgene_id=as.character(gene_id_ms)), by="entrezgene_id")