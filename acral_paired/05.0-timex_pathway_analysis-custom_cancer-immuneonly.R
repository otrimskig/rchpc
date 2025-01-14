source("libs.R")

#cudc rna seq timex analysis 
#timex does not care about comparisons, so any group loaded is fine 
#needs RPKM data in a format where the rownames are the gene names (gene symbol) and the rpkms are a matrix 


library(tidyverse)
library(dtplyr)
library(GSVA)


#read in data to be used. format such that is in a wide data matrix,
#with rpkms as values, sample names as column titles,
#and rownames as unique(!!), human(!!) gene names. 
rpkms<-readRDS("acral_paired/ds/v02-filtered_rpkms.rds")
name_map<-readRDS("acral_paired/ds/name_map.rds")
gene_stats<-readRDS("acral_paired/ds/gene_stats.rds")



#need to get unique human gene names. 
#must be unique otherwise won't be able to convert to rownames. 
#where duplicates exist, choose the one to use based
#on which had the higher mean_rpkm. Use that to select the 
#corresponding gene_id_ms, which should be unique for each. 
#can then use that to filter going forward.
ms_gene_ids_unique_hu<-gene_stats%>%
  filter(!is.na(gene_name_hu))%>%
  group_by(gene_name_hu)%>%
  arrange(desc(mean_rpkm))%>%
  slice(1)%>%ungroup()%>%
  select(gene_id_ms)


rpkms2<-rpkms%>%
  semi_join(ms_gene_ids_unique_hu)%>%
  left_join(name_map%>%select(sample_id, mouse_num_type_clean))%>%
  select(-sample_id)%>%
  select(gene_name_hu, rpkm, mouse_num_type_clean)%>%
  pivot_wider(values_from = rpkm, names_from = mouse_num_type_clean)%>%
  column_to_rownames("gene_name_hu")


mat<-rpkms2%>%data.matrix()


#get, as a named list, all the gene signatures.
#they are named list of pathway genes.

load("timex/ds/allSignatures.rda")  

names(Hallmark) <-sub("^HALLMARK_", "hallmark_", names(Hallmark))
names(TIMEx) <- paste0("timex_", names(TIMEx))
names(Immune_sig) <- paste0("immune_", names(Immune_sig))

onco<-qusage::read.gmt("timex/ds/c6.all.v2024.1.Hs.symbols.gmt")
names(onco) <- paste0("onco_", names(onco))


#create signatures list from elements you want to include in analysis.
Signature_list <- c(Hallmark, kegg, TIMEx, Immune_sig, onco)

#make parameters list.
ssgeaP<-ssgseaParam(mat, Signature_list)

# ssgeaP<-gsvaParam(mat, Signature_list)




#normalize gene set expression. suppress warnings to ignore
#sets that are only 1 gene. 


gng_gsva<- suppressWarnings(gsva(ssgeaP))
gng_ssgsea_u <- t(t(gng_gsva))


#sometimes strips a few signatures from the analysis not sure why. this reveals what is unique 
namesGSEA<-as.list(rownames(gng_ssgsea_u))
gseavector<-unlist(namesGSEA)
signaturevector<-unlist(names(Signature_list))





gng_ssgsea_z <- t(scale(t(gng_gsva)))

#sometimes strips a few signatures from the analysis not sure why. this reveals what is unique 
namesGSEA<-as.list(rownames(gng_ssgsea_z))

gseavector<-unlist(namesGSEA)


signaturevector<-unlist(names(Signature_list))

signaturevector<- stringr::str_remove_all(signaturevector, "HALLMARK_")


notincluded<-setdiff(signaturevector, gseavector)%>%
  tibble(signature_name=., included="0")




stop("check file save locations")
# 
saveRDS(gng_ssgsea_u, "acral_paired/ds/gsva_u.rds")
saveRDS(gng_ssgsea_z, "acral_paired/ds/gsva_z.rds")




tibble(signature_name=signaturevector, included="1")%>%
  full_join(notincluded)%>%
  rename(gsva_name=signature_name)%>%
  mutate(list_name = str_split(gsva_name, "_", simplify = TRUE)[, 1])%>%
  mutate(pathway_name = map2_chr(list_name, gsva_name, ~sub(paste0("^", .x, "_"), "", .y)))%>%
  mutate(pathway_name_clean=toupper(pathway_name))%>%
  mutate(pathway_name_clean=gsub("_", " ", pathway_name_clean))%>%

  
  
  saveRDS("acral_paired/ds/gsva_pathway_names.rds")
  
  
