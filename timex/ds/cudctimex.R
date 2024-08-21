
#cudc rna seq timex analysis 
# timex does not care about comparisons, so any group loaded is fine 
#needs RPKM data in a format where the rownames are the gene names (gene symbol) and the rpkms are a matrix 
library(tidyverse)
library(GSVA)
library(ComplexHeatmap)
library(colorRamp2)

edgeRout<- readr::read_tsv(file = "C:/Users/Michael/OneDrive - University of Utah/Hawse/CUDC101RNAseq/REDO/eDGEevVec.txt")


RPKMcol<-edgeRout[,25:36]
# cleanup to make RPKMcol look normal 
RPKMcolclean<-cbind(RPKMcol[,10:12], RPKMcol[,4:6], RPKMcol[,7:9], RPKMcol[,1:3])
#add the gene symbols to the first column 
cudctimex<-cbind(edgeRout$external_gene_name, RPKMcolclean)
cudctimex2<-cudctimex %>%
  distinct(`edgeRout$external_gene_name`, .keep_all = TRUE) %>%
  na.omit() 
  
rownames(cudctimex2)<-cudctimex2$`edgeRout$external_gene_name`
cudctimex2<-cudctimex2[-1]

toMat2 <- function(input){
  genes <- input$`edgeRout$external_gene_name`
  input <- input[, -1]
  input <- as.matrix(input)
  rownames(input) <- genes
  input
}


gng_mat <- toMat2(cudctimex2)

load("C:/Users/Michael/OneDrive - University of Utah/Holmen/projects/timex/allSignatures.rda")

Signature_list <- c(Hallmark, kegg, TIMEx, Immune_sig)
#as.matrix() works 
gng_ssgsea <- GSVA::gsva(as.matrix(cudctimex2), gset.idx.list = Signature_list, method = "ssgsea")
rownames(gng_ssgsea) <- stringr::str_remove_all(rownames(gng_ssgsea), "HALLMARK_")

gng_ssgsea_z <- t(scale(t(gng_ssgsea)))

#sometimes strips a few signatures from the analysis not sure why. this reveals what is unique 
namesGSEA<-as.list(rownames(gng_ssgsea_z))
gseavector<-unlist(namesGSEA)
signaturevector<-unlist(names(Signature_list))
signaturevector<- stringr::str_remove_all(signaturevector, "HALLMARK_")
notincluded<-setdiff(signaturevector, gseavector)


sig_tracker_df <- data.frame("Signature" = gseavector,
                             "Set" = c( 
                               rep("Hallmark", times = 50), 
                               rep("KEGG", times = 186),
                               rep("TIMEx", times = 37),
                               rep("Immune", times = 46))) #48 originally updated with setdiff 




gng_heat <- Heatmap(gng_ssgsea_z,
                    name = "cudcRNAseq",
                    #column_split = factor(genepuree2$Growth, levels = c("Grow", "nogrow")),
                    cluster_columns= FALSE,
                    cluster_column_slices = FALSE,
                    #top_annotation=genepuree3,
                    row_split = factor(sig_tracker_df$Set, levels = c(
                      "Hallmark",
                      "KEGG",
                      "TIMEx",
                      "Immune")), #sometimes fails not sure why 
                    cluster_row_slices = FALSE,
                    row_names_gp = gpar(fontsize = 6),
                    column_names_gp = gpar(fontsize =2),
                    column_title_gp = gpar(font = 2, fontsize = 60),
                    row_title_gp = gpar(font = 2, fontsize = 40),
                    heatmap_legend_param = list(
                      labels_gp = gpar(fontsize = 20),
                      title_gp = gpar(fontsize = 24)
                    ))
plot(gng_heat)
# Save PDF 30 x 60 inches #







