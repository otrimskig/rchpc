
library(tidyverse)
library(GSVA)
library(ComplexHeatmap)
library(colorRamp2)

gng <- readr::read_tsv("TPM_dod_grownogrow_for_cibersortx.txt")

toMat2 <- function(input){
  genes <- input$Gene
  input <- input[, -1]
  input <- as.matrix(input)
  rownames(input) <- genes
  input
}

gng_mat <- toMat2(gng)

load("allSignatures.rda")

Signature_list <- c(Hallmark, kegg, TIMEx, Immune_sig)

gng_ssgsea <- GSVA::gsva(gng_mat, gset.idx.list = Signature_list, method = "ssgsea")
rownames(gng_ssgsea) <- stringr::str_remove_all(rownames(gng_ssgsea), "HALLMARK_")

gng_ssgsea_z <- t(scale(t(gng_ssgsea)))

sig_tracker_df <- data.frame("Signature" = names(Signature_list),
                             "Set" = c( 
                               rep("Hallmark", times = 50), 
                               rep("KEGG", times = 186),
                               rep("TIMEx", times = 37),
                               rep("Immune", times = 48)))

#genepuree3<-HeatmapAnnotation(df=genepuree2, show_legend=TRUE)

#sub_samp_ordered<-gng_ssgsea_z[,rownames(genepuree2)]

gng_heat <- Heatmap(gng_ssgsea_z,
                    name = "DOD GrownoGrow",
                    #column_split = factor(genepuree2$Growth, levels = c("Grow", "nogrow")),
                    cluster_columns= FALSE,
                    cluster_column_slices = FALSE,
                    top_annotation=genepuree3,
                    row_split = factor(sig_tracker_df$Set, levels = c(
                      "Hallmark",
                      "KEGG",
                      "TIMEx",
                      "Immune")),
                    cluster_row_slices = FALSE,
                    row_names_gp = gpar(fontsize = 6),
                    column_names_gp = gpar(fontsize =8),
                    column_title_gp = gpar(font = 2, fontsize = 60),
                    row_title_gp = gpar(font = 2, fontsize = 40),
                    heatmap_legend_param = list(
                      labels_gp = gpar(fontsize = 20),
                      title_gp = gpar(fontsize = 24)
                    ))

plot(gng_heat)
# Save PDF 30 x 60 inches #

