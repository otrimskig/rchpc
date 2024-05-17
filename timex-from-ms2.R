# if (!exists("n.cores")) {
#   
#   "initilizing cores..."
#   n.cores <- parallel::detectCores() - 1
#   my.cluster <- parallel::makeCluster(
#     n.cores, 
#     type = "PSOCK"
#   )
#   doParallel::registerDoParallel(cl = my.cluster)
#   
#   #check if it is registered (optional)
#   foreach::getDoParRegistered()
#   
#   "parallel cores initialized."
#   
# }
# 
# 





source("libs.R")

#cudc rna seq timex analysis 
# timex does not care about comparisons, so any group loaded is fine 
#needs RPKM data in a format where the rownames are the gene names (gene symbol) and the rpkms are a matrix 
library(tidyverse)
library(dtplyr)
suppressPackageStartupMessages(library(ComplexHeatmap))


##########################################
# all<-readRDS("ds/v08-rpkms_w_meta_ms_hu.rds")%>%
#   filter(!is.na(gene_name_ms))%>%
#   filter(!is.na(gene_name_hu))%>%
# 
#   select(gene_name_hu, mouse_num, rpkm)%>%
#   mutate(mouse_num=paste0("x", mouse_num))%>%
#   
#   group_by(gene_name_hu, mouse_num)%>%slice(1)%>%ungroup()
# 
#   
# 
# all2<-all%>%pivot_wider(names_from = mouse_num, values_from = rpkm)
#   
#   
# saveRDS(all2, "ds/v09-hu-gene-rpkms-wide.rds")

#####

mat<-readRDS("ds/v09-hu-gene-rpkms-wide.rds")%>%
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
# saveRDS(gng_ssgsea_u, "ds/gsva_u.rds")
# saveRDS(gng_ssgsea_z, "ds/gsva_z.rds")
# 
# saveRDS(sig_tracker_df_u, "ds/gsva_sig_u.rds")
# saveRDS(sig_tracker_df_z, "ds/gsva_sig_z.rds")
# 
# stop()
gng_heat <- Heatmap(gng_ssgsea_z,
                    #name = "cudcRNAseq",
                    #column_split = factor(genepuree2$Growth, levels = c("Grow", "nogrow")),
                    #cluster_columns= FALSE,
                    cluster_column_slices = FALSE,
                    #top_annotation=genepuree3,
                    row_split = factor(sig_tracker_df$Set, levels = c(
                      "Hallmark",
                      "KEGG",
                      "TIMEx",
                      "Immune")), #sometimes fails not sure why 
                    #cluster_row_slices = FALSE,
                    #row_names_gp = gpar(fontsize = 6),
                    #column_names_gp = gpar(fontsize =2),
                    column_title_gp = gpar(font = 2, fontsize = 60),
                    row_title_gp = gpar(font = 2, fontsize = 40),
                    show_heatmap_legend = FALSE,
                      
                    
                   width = unit(15, "in"),  
                    heatmap_height = unit(65, "in")
                    
                  )
gh<-grid.grabExpr(draw(gng_heat))


ggsave("plots/hm-timex-all-rpkm-unscaled.pdf",
       plot=gh,
       
       scale = 1,
       dpi=600,
       width = 28,
       height = 68,
       unit="in",
       limitsize = FALSE
       
)

