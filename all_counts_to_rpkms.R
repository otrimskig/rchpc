library(tidyverse)
library(dtplyr)
library(biomaRt)

df<-readRDS("ds/mv-00-all_counts.rds")%>%
  rename(gene_id_ms=GeneID)%>%
  rename(gene_len_ms=Length)%>%
  rename_at(vars(contains("merged")), ~ sub("_merged.bam", "", .))






human <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

ms_gene_info <- getBM(attributes = c("entrezgene_id", "external_gene_name"),
                      filters = "entrezgene_id",
                      values = df$gene_id_ms,
                      mart = mouse) %>%
  rename(gene_id_ms = entrezgene_id, gene_name_ms = external_gene_name)


subset_genes <- trimws(ms_gene_info$gene_id_ms[1:10])


# ms_gene_info<-getBM(attributes = c("entrezgene_id", "external_gene_name"),
#                  filters="entrezgene_id",
#                  values = df$gene_id_ms, 
#                  mart = mouse)%>%
#   rename(gene_id_ms=entrezgene_id, gene_name_ms=external_gene_name)

minimal_genes <- c("BRCA1", "TP53", "EGFR")
hu_gene_info <- getLDS(attributes = c("hgnc_symbol"),
                       filters = "hgnc_symbol",
                       values = minimal_genes,
                       mart = human,
                       attributesL = c("mgi_symbol"),
                       martL = mouse,
                       uniqueRows = TRUE)

hu_gene_info <- getLDS(attributes = c("hgnc_symbol"),
                       filters = "hgnc_symbol",
                       values = subset_genes,
                       mart = human,
                       
                       attributesL = c("entrezgene_id"),
                       martL = mouse,
                       uniqueRows = TRUE)


# hu_gene_info<- getLDS(attributes = c("hgnc_symbol"), 
#                    filters = "hgnc_symbol", 
#                    values = ms_gene_info$gene_name_ms, 
#                    mart = human, 
#                    attributesL = c("mgi_symbol"), 
#                    martL = mouse, uniqueRows=T)










# gene_info<-gene_info%>%
#   rename(gene_id_ms=entrezgene_id, gene_name_ms=external_gene_name)
#   
# gene_info2<-gene_info%>%group_by(gene_id_ms)%>%slice(1)%>%ungroup()%>%
#   group_by(gene_name_ms)%>%slice(1)%>%ungroup()



# name_counts<-gene_info2%>%
#   group_by(gene_name_ms)%>%
#   summarise(unique_gene_names=n())
# 
# 
# 
# 
# 
# 
# df3<-df2%>%
#   left_join(gene_info2)






stop("done")







df4<-df3%>%
  pivot_longer(cols=contains("X23908"), names_to="sample_id", values_to="read_count")%>%
  group_by(sample_id)%>%
 
  mutate(rpkm = read_count/(gene_len_ms/1000 * (sum(read_count)/1E6)))




head(df4)%>%view()
