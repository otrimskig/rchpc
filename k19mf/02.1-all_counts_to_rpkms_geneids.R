library(tidyverse)
library(dtplyr)
library(biomaRt)

#relative to rchpc proj
proj_dir<-"k19mf/"


df<-readRDS(paste0(proj_dir, "ds/vm-00-all_counts_sample_names_fixed.rds"))%>%
  rename(gene_len_ms=length)%>%
  rename(gene_id_ms=gene_id)
  #rename_at(vars(contains("merged")), ~ sub("_merged.bam", "", .))






human <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org/")
mouse <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host="https://dec2021.archive.ensembl.org/")

ms_gene_info <- getBM(attributes = c("entrezgene_id", "external_gene_name"),
                      filters = "entrezgene_id",
                      values = df$gene_id_ms,
                      mart = mouse) %>%
  rename(gene_id_ms = entrezgene_id, gene_name_ms = external_gene_name)



human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

hu_gene_info = getLDS(attributes = c("hgnc_symbol"), 
                      filters = "hgnc_symbol", 
                      values = ms_gene_info$gene_name_ms , 
                      mart = human, attributesL = c("mgi_symbol"), 
                      martL = mouse, uniqueRows=T)%>%
  rename(gene_name_ms=MGI.symbol, gene_name_hu= HGNC.symbol)




gene_info<-ms_gene_info%>%
  left_join(hu_gene_info)



unique_genes<-gene_info%>%
  mutate(gene_name_ms=trimws(gene_name_ms))%>%
  mutate(gene_name_ms=na_if(gene_name_ms, ""))%>%
  group_by(gene_id_ms)%>%
  arrange(gene_name_hu, gene_name_ms)%>%
  slice(1)%>%
  ungroup()





df2<-df%>%
  left_join(unique_genes)%>%
  relocate(gene_id_ms, gene_len_ms, gene_name_ms, gene_name_hu)%>%
  
  # pivot_longer(cols=contains("X23908"), names_to="sample_id", values_to="read_count")%>%
 
   pivot_longer(cols=5:last_col(), names_to="sample_id", values_to="read_count")%>%
  
  group_by(sample_id)%>%
 
  mutate(rpkm = read_count/(gene_len_ms/1000 * (sum(read_count)/1E6)))%>%
  group_by(gene_id_ms)
  




saveRDS(df2, paste0(proj_dir, "ds/vm-01-gene_id_hu_m_rpkms.rds"))
