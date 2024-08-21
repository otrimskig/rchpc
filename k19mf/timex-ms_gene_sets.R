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

onco<-qusage::read.gmt("timex/ds/c6.all.v2024.1.Hs.symbols.gmt")
names(onco) <- paste0("onco_", names(onco))


mhc2<-Signature_list$MHC_II


exp_data<-readRDS("k19mf/ds/vm-h-01-rpkms_wide_human.rds")%>%
  filter(gene_name_hu %in% mhc2)



exp_data2<-readRDS("k19mf/ds/gene_stats.rds")%>%
  filter(gene_name_hu %in% mhc2)





#create signatures list from elements you want to include in analysis.
Signature_list <- c(Hallmark, kegg, TIMEx, Immune_sig, onco)

#make parameters list.
ssgeaP<-ssgseaParam(mat, Signature_list)

gsvaP<-gsvaParam(mat, Signature_list)

#normalize gene set expression. suppress warnings to ignore
#sets that are only 1 gene. 




gng_ssgsea<- gsva(ssgeaP)

gng_ssgsea<- gsva(gsvaP)




suppressWarnings(gng_ssgsea<- gsva(ssgeaP))

###
rownames(gng_ssgsea)


rownames(gng_ssgsea) <- stringr::str_remove_all(rownames(gng_ssgsea), "HALLMARK_")
gng_ssgsea_u <- t(t(gng_ssgsea))


#sometimes strips a few signatures from the analysis not sure why. this reveals what is unique 
namesGSEA<-as.list(rownames(gng_ssgsea_u))
gseavector<-unlist(namesGSEA)


gseavector




signaturevector<-unlist(names(Signature_list))

signaturevector

signaturevector<- stringr::str_remove_all(signaturevector, "HALLMARK_")




notincluded<-setdiff(signaturevector, gseavector)

notincluded






sig_tracker_df_u <- data.frame("Signature" = gseavector,
                               "Set" = c( 
                                 rep("Hallmark", times = 50), 
                                 rep("KEGG", times = 186),
                                 rep("TIMEx", times = 37),
                                 rep("Immune", times = 46),
                                 rep("onco", times = 189)
                                 )
                               ) #48 originally updated with setdiff 


gng_ssgsea_z <- t(scale(t(gng_ssgsea)))

#sometimes strips a few signatures from the analysis not sure why. this reveals what is unique 
namesGSEA<-as.list(rownames(gng_ssgsea_z))
gseavector<-unlist(namesGSEA)
signaturevector<-unlist(names(Signature_list))
signaturevector<- stringr::str_remove_all(signaturevector, "HALLMARK_")
notincluded<-setdiff(signaturevector, gseavector)



sig_tracker_df_z <- data.frame("Signature" = gseavector,
                               "Set" = c( 
                                 rep("Hallmark", times = 50), 
                                 rep("KEGG", times = 186),
                                 rep("TIMEx", times = 37),
                                 rep("Immune", times = 46),
                                 rep("onco", times = 189)
                                )
                               ) #48 originally updated with setdiff 



stop("check file save locations")
# 
saveRDS(gng_ssgsea_u, "k19mf/ds/gsva_u-onco.rds")
saveRDS(gng_ssgsea_z, "k19mf/ds/gsva_z-onco.rds")
# 
saveRDS(sig_tracker_df_u, "k19mf/ds/gsva_sig_u-onco.rds")
saveRDS(sig_tracker_df_z, "k19mf/ds/gsva_sig_z-onco.rds")










