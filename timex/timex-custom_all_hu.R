source("libs.R")

#cudc rna seq timex analysis 
#timex does not care about comparisons, so any group loaded is fine 
#needs RPKM data in a format where the rownames are the gene names (gene symbol) and the rpkms are a matrix 


library(tidyverse)
library(dtplyr)
library(GSVA)



mat<-readRDS("nf1g/ds/vm-h-01-rpkms_wide_human.rds")%>%
  column_to_rownames("gene_name_hu")%>%
  data.matrix()


# 
# load("timex/ds/allSignatures.rda")  

onco<-qusage::read.gmt("timex/ds/c6.all.v2024.1.Hs.symbols.gmt")
names(onco) <- paste0("onco_", names(onco))






Signature_list <- c(Hallmark, kegg, TIMEx, Immune_sig, onco)




ssgeaP<-ssgseaParam(mat, Signature_list)




suppressWarnings(gng_ssgsea<- gsva(ssgeaP))
rownames(gng_ssgsea) <- stringr::str_remove_all(rownames(gng_ssgsea), "HALLMARK_")
gng_ssgsea_u <- t(t(gng_ssgsea))


#sometimes strips a few signatures from the analysis not sure why. this reveals what is unique 
namesGSEA<-as.list(rownames(gng_ssgsea_u))
gseavector<-unlist(namesGSEA)
signaturevector<-unlist(names(Signature_list))
signaturevector<- stringr::str_remove_all(signaturevector, "HALLMARK_")
notincluded<-setdiff(signaturevector, gseavector)


sig_tracker_df_u <- data.frame("Signature" = gseavector,
                               "Set" = c( 
                                 rep("Hallmark", times = 50), 
                                 rep("KEGG", times = 186),
                                 rep("TIMEx", times = 37),
                                 rep("Immune", times = 46),
                                 rep("onco", times = 189),
                                 rep("gsea_all", times=34837))) #48 originally updated with setdiff 


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
                                 rep("onco", times = 189),
                                 rep("gsea_all", times=34837))) #48 originally updated with setdiff 


# 
saveRDS(gng_ssgsea_u, "nf1g/ds/gsva_u-onco.rds")
saveRDS(gng_ssgsea_z, "nf1g/ds/gsva_z-onco.rds")
# 
saveRDS(sig_tracker_df_u, "nf1g/ds/gsva_sig_u-onco.rds")
saveRDS(sig_tracker_df_z, "nf1g/ds/gsva_sig_z-onco.rds")










