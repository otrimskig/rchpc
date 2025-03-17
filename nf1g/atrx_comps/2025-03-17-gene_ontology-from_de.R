source("libs.R")
library(tidyverse)
library(dtplyr)
library(ggplot2)
library(biomaRt)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)


#get vector of all .rds files containing dexp analyses.
all_dexps<-list.files("nf1g/dexps/atrx_comps", full.names = T, pattern="^dexp.*\\.rds$")







gene_stats<-readRDS("nf1g/ds/gene_stats.rds")
gene_list<-gene_stats%>%dplyr::select(gene_name_ms)%>%drop_na()%>%
  unique()%>%pull()

df0<-readRDS(all_dexps[1])



# Connect to Ensembl Mouse database
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")



# Retrieve Gene Ontology (GO) terms
go_info <- getBM(
  attributes = c("external_gene_name", "go_id", "name_1006"),  # GO ID and description
  filters = "external_gene_name",
  values = gene_list,
  mart = mart
)
go_info<-go_info%>%
  as_tibble()%>%
  dplyr::rename(Gene = "external_gene_name", GO_ID = "go_id", GO_Description = "name_1006")

# Retrieve Reactome Pathways (Alternative to KEGG)
reactome_info <- getBM(
  attributes = c("external_gene_name", "reactome"),
  filters = "external_gene_name",
  values = gene_list,
  mart = mart
) %>%
  dplyr::rename(Gene = external_gene_name, Reactome_Pathway = reactome)

# Convert Gene Symbols to Entrez IDs for KEGG analysis
entrez_ids <- mapIds(
  org.Mm.eg.db,
  keys = gene_list,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
) %>% na.omit()

# Perform KEGG Pathway Enrichment using clusterProfiler
kegg_results <- enrichKEGG(
  gene = entrez_ids,
  organism = "mmu",  # Mouse KEGG prefix
  pvalueCutoff = 0.05
)

# Convert KEGG results to a tibble
kegg_info <- as_tibble(kegg_results@result) %>%
  dplyr::select(ID, Description, geneID) %>%
  dplyr::rename(KEGG_Pathway_ID = ID, KEGG_Description = Description, Genes_Involved = geneID)

# Retrieve Pfam (Protein Family) annotations
pfam_info <- getBM(
  attributes = c("external_gene_name", "pfam"),
  filters = "external_gene_name",
  values = gene_list,
  mart = mart
) %>%
  dplyr::rename(Gene = external_gene_name, Pfam_ID = pfam)

# Retrieve InterPro Domains (Protein Domains)
interpro_info <- getBM(
  attributes = c("external_gene_name", "interpro", "interpro_description"),
  filters = "external_gene_name",
  values = gene_list,
  mart = mart
) %>%
  dplyr::rename(Gene = external_gene_name, InterPro_ID = interpro, InterPro_Description = interpro_description)

# Ensure each dataset has a "Gene" column and replace empty entries with NA
combined_info <- list(
  GO = go_info,
  Reactome = reactome_info,
  KEGG = kegg_info,
  Pfam = pfam_info,
  InterPro = interpro_info
)



combined_info2<-lapply(combined_info, as_tibble)


# Expand the KEGG dataset so that each gene is in its own row
KEGG_expanded <- combined_info2$KEGG %>%
  separate_rows(Genes_Involved, sep = "/") %>%  # Split genes by "/"
  mutate(Genes_Involved = trimws(Genes_Involved))%>%
  dplyr::rename(ensembl_id = Genes_Involved)# Remove extra spaces around genes

# Now, KEGG_expanded will have one row per gene for each pathway
# You can join this back into the combined_info list by "Gene"

# Replace the original KEGG dataset with the expanded one in the combined_info list
combined_info2$KEGG <- KEGG_expanded


combined_info2$KEGG <- combined_info2$KEGG%>%
  dplyr::rename(gene_id_ms=ensembl_id)%>%
  left_join(gene_stats%>%mutate(gene_id_ms=as.character(gene_id_ms)))%>%
  select(starts_with("KEGG"), gene_name_ms)%>%
  filter(!is.na(gene_name_ms))%>%
  dplyr::rename(Gene=gene_name_ms)





combined_info3<-combined_info2%>%
  
  map(~.x %>%
        # Rename if needed
        mutate_all(~na_if(., "") %>% na_if(" "))) 






# Now you can try to join the datasets as before
combined_info4 <- lapply(combined_info3, as_tibble)



combined_info5<-gene_stats%>%
  select(gene_name_ms)%>%
  unique()%>%
  drop_na()%>%
  dplyr::rename(Gene=gene_name_ms)%>%
  left_join(combined_info4$GO)%>%
  drop_na()%>%
  left_join(combined_info4$Reactome)%>%
  drop_na()%>%
  left_join(combined_info4$KEGG)%>%
  drop_na()%>%
  left_join(combined_info4$)





saveRDS(combined_info3, "nf1g/ds/gene_combined_info_go.rds")

# Print combined results
view(sample_n(combined_info4, 1000))
