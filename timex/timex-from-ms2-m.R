source("libs.R")

#cudc rna seq timex analysis 
# timex does not care about comparisons, so any group loaded is fine 
#needs RPKM data in a format where the rownames are the gene names (gene symbol) and the rpkms are a matrix 
library(tidyverse)
library(dtplyr)
suppressPackageStartupMessages(library(ComplexHeatmap))

#####

data<-readRDS("ds/vm-02-filtered_rpkms.rds")

unique_hu_genes<-readRDS("ds/gene_stats.rds")%>%
  filter(!is.na(gene_name_hu))%>%
  group_by(gene_name_hu)%>%
  arrange(desc(mean_rpkm))%>%
  slice(1)%>%ungroup()


mouse_nums<-readRDS("ds/v10-per_sample_updated.rds")%>%
  select(mouse_num, sample_id)


data_hu_only<-data%>%
  semi_join(gene_stats, by="gene_id_ms")%>%
  left_join(mouse_nums)%>%
  select(-sample_id)%>%
  select(gene_name_hu, rpkm, mouse_num)%>%
  pivot_wider(names_from = mouse_num, values_from=rpkm)


saveRDS(data_hu_only, "ds/vm-h-01-rpkms_wide_human.rds")

mat<-data_hu_only%>%
  column_to_rownames("gene_name_hu")%>%
  data.matrix()
  





library(GSVA)

load("timex/allSignatures.rda")  
Signature_list <- c(Hallmark, kegg, TIMEx, Immune_sig)



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
                               rep("Immune", times = 46))) #48 originally updated with setdiff 


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
                                 rep("Immune", times = 46))) #48 originally updated with setdiff 


# 
saveRDS(gng_ssgsea_u, "ds/gsva_u.rds")
saveRDS(gng_ssgsea_z, "ds/gsva_z.rds")
# 
saveRDS(sig_tracker_df_u, "ds/gsva_sig_u.rds")
saveRDS(sig_tracker_df_z, "ds/gsva_sig_z.rds")


