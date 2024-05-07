library(tidyverse)
library(biomaRt)


ensembl <- useEnsembl(biomart = "genes")

ensembl<-useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)




gene_counts<-read.table("23908R/v01-all_counts.txt", header = TRUE)


info<-c("entrezgene_id", "external_gene_name")



gene_info<-getBM(attributes = info,
      filter="entrezgene_id",
      values = gene_counts$GeneID, 
      mart = ensembl)


gene_info%>%
  count(external_gene_name)%>%
  view()



name_counts<-gene_info%>%
  group_by(entrezgene_id)%>%
  summarise(unique_gene_names=n())


unique<-name_counts%>%
  left_join(gene_info)%>%
  group_by(entrezgene_id)%>%
  slice(1)%>%
  ungroup()
  





fg<-gene_counts%>%
  left_join(unique, by=c("GeneID"="entrezgene_id"))%>%
  relocate(unique_gene_names, external_gene_name, .after="Length")%>%
  janitor::clean_names()





clean_column_names <- function(col_name) {
  str_extract(col_name, "^[^_]+")
}

# Get column names starting from the 5th column
cols_to_clean <- names(fg)[5:length(fg)]

# Clean column names
cleaned_names <- cols_to_clean %>%
  map(~ clean_column_names(.))

# Rename the columns
names(fg)[5:length(fg)] <- cleaned_names



fg2<-fg%>%
  rename(x23908x01=x23908x1,
         x23908x02=x23908x2,
         x23908x03=x23908x3,
         x23908x04=x23908x4,
         x23908x05=x23908x5,
         x23908x06=x23908x6,
         x23908x07=x23908x7,
         x23908x08=x23908x8,
         x23908x09=x23908x9)



fg2<-fg2[, order(names(fg2))]







saveRDS(fg2, "23908R/v02-counts_with_geneIDS.rds")




