source("libs.R")

#cudc rna seq timex analysis 
#timex does not care about comparisons, so any group loaded is fine 
#needs RPKM data in a format where the rownames are the gene names (gene symbol) and the rpkms are a matrix 


library(tidyverse)
library(dtplyr)
library(GSVA)


readRDS("nf1g/ds/v10-per_sample_updated.rds")

ms_mat<-readRDS("nf1g/ds/vm-02-filtered_rpkms.rds")%>%
  select(gene_id_ms, sample_id, rpkm)%>%
  pivot_wider(names_from = sample_id, values_from = rpkm)%>%
  column_to_rownames("gene_id_ms")%>%
  data.matrix()

ms_pathways<-qusage::read.gmt("timex/ds/m8.all.v2024.1.Mm.entrez-celltypesig.gmt")

gng_ssgsea<-gsva(ssgseaParam(ms_mat, ms_pathways))%>%
  suppressWarnings()



gng_ssgsea_u <- t(t(gng_ssgsea))
gng_ssgsea_z <- t(scale(t(gng_ssgsea)))


namesGSEA<-as.list(rownames(gng_ssgsea_u))
gseavector<-unlist(namesGSEA)


signaturevector<-unlist(names(ms_pathways))


notincluded<-setdiff(signaturevector, gseavector)



saveRDS(gng_ssgsea_u, "nf1g/ds/gsva_u_celltypes.rds")
saveRDS(gng_ssgsea_z, "nf1g/ds/gsva_z_celltypes.rds")
