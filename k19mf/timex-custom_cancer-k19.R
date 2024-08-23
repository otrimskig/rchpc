source("libs.R")

#cudc rna seq timex analysis 
#timex does not care about comparisons, so any group loaded is fine 
#needs RPKM data in a format where the rownames are the gene names (gene symbol) and the rpkms are a matrix 


library(tidyverse)
library(dtplyr)
library(GSVA)


#read in data to be used. format such that is in a wide data matrix,
#with rpkms as values, sample names as column titles,
#and rownames as unique, human(!) gene names. 
mat<-readRDS("k19mf/ds/vm-h-01-rpkms_wide_human.rds")%>%
  column_to_rownames("gene_name_hu")%>%
  data.matrix()



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
saveRDS(gng_ssgsea_u, "k19mf/ds/gsva_u-onco.rds")
saveRDS(gng_ssgsea_z, "k19mf/ds/gsva_z-onco.rds")




tibble(signature_name=signaturevector, included="1")%>%
  full_join(notincluded)%>%
  rename(gsva_name=signature_name)%>%
  mutate(list_name = str_split(gsva_name, "_", simplify = TRUE)[, 1])%>%
  mutate(pathway_name = map2_chr(list_name, gsva_name, ~sub(paste0("^", .x, "_"), "", .y)))%>%
  mutate(pathway_name_clean=toupper(pathway_name))%>%
  mutate(pathway_name_clean=gsub("_", " ", pathway_name_clean))%>%

  
  
  saveRDS("k19mf/ds/gsva_pathway_names.rds")
  
  
