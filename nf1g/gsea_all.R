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




gsea_all_sets<-qusage::read.gmt("timex/ds/msigdb.v2024.1.Hs.symbols.gmt")




Signature_list <- c(gsea_all_sets)




ssgeaP<-ssgseaParam(mat, Signature_list)




suppressWarnings(gng_ssgsea<- gsva(ssgeaP))
rownames(gng_ssgsea) <- stringr::str_remove_all(rownames(gng_ssgsea), "HALLMARK_")
gng_ssgsea_u <- t(t(gng_ssgsea))


#sometimes strips a few signatures from the analysis not sure why. this reveals what is unique 
namesGSEA<-as.list(rownames(gng_ssgsea_u))
gseavector<-unlist(namesGSEA)
signaturevector<-unlist(names(Signature_list))

notincluded<-setdiff(signaturevector, gseavector)




gng_ssgsea_z <- t(scale(t(gng_ssgsea)))

#sometimes strips a few signatures from the analysis not sure why. this reveals what is unique 
namesGSEA<-as.list(rownames(gng_ssgsea_z))
gseavector<-unlist(namesGSEA)
signaturevector<-unlist(names(Signature_list))

notincluded<-setdiff(signaturevector, gseavector)






# 
saveRDS(gng_ssgsea_u, "nf1g/ds/gsva_all_u.rds")
saveRDS(gng_ssgsea_z, "nf1g/ds/gsva_all_z.rds")
# 













