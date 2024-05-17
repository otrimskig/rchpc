#combination heatmap with "z-score" heat. 

source("libs.R")
library(tidyverse)
library(dtplyr)
suppressPackageStartupMessages(library(ComplexHeatmap))
library(gridtext)


#load relevant created from gsva analysis. 
gsva_u_sig<-readRDS("ds/gsva_sig_u.rds")
gsva_z_sig<-readRDS("ds/gsva_sig_z.rds")
gsva_u<-readRDS("ds/gsva_u.rds")
gsva_z<-readRDS("ds/gsva_z.rds")







hm_z <- Heatmap(gsva_z,
                
                
                column_title = gt_render(
                  paste0("<span style='font-size:25pt'>Expression Levels Normalized per gene <br>(Relative expression)</span>"), 
                  r = unit(2, "pt")),
                    #name = "cudcRNAseq",
                    #column_split = factor(genepuree2$Growth, levels = c("Grow", "nogrow")),
                    #cluster_columns= FALSE,
                    cluster_column_slices = FALSE,
                    #top_annotation=genepuree3,
                    row_split = factor(gsva_z_sig$Set, levels = c(
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
                    
                 
                    width = unit(12, "in"),
                    heatmap_height = unit(65, "in")

)
hm_u <- Heatmap(gsva_u,
                column_title = gt_render(
                  paste0("<span style='font-size:25pt'>Absolute <br>expression</span>"), 
                  r = unit(2, "pt")),
                #name = "cudcRNAseq",
                #column_split = factor(genepuree2$Growth, levels = c("Grow", "nogrow")),
                cluster_columns= FALSE,
                column_order=column_order(hm_z),
                show_column_dend = TRUE,
                cluster_column_slices = FALSE,
                #top_annotation=genepuree3,
                row_split = factor(gsva_u_sig$Set, levels = c(
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
                
                
                width = unit(2, "in"),  
                heatmap_height = unit(65, "in")
                
)



combo<-hm_u+hm_z



gh<-grid.grabExpr(draw(combo))


ggsave("plots/hm-timex-all-rpkm-combo-clustered_abs-local.pdf",
       plot=gh,
       
       scale = 1,
       dpi=600,
       width = 45,
       height = 100,
       unit="in",
       limitsize = FALSE
       
)


#######################################################################



hm_z2 <- Heatmap(gsva_z,
                
                
                column_title = gt_render(
                  paste0("<span style='font-size:25pt'>Expression Levels Normalized per gene <br>(Relative expression)</span>"), 
                  r = unit(2, "pt")),
                #name = "cudcRNAseq",
                #column_split = factor(genepuree2$Growth, levels = c("Grow", "nogrow")),
                #cluster_columns= FALSE,
                cluster_column_slices = FALSE,
                #top_annotation=genepuree3,
                row_split = factor(gsva_z_sig$Set, levels = c(
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
                
                
                width = unit(12, "in"),
                heatmap_height = unit(65, "in")
                
)

row_order(hm_z2)

hm_u2 <- Heatmap(gsva_u,
                column_title = gt_render(
                  paste0("<span style='font-size:25pt'>Absolute <br>expression</span>"), 
                  r = unit(2, "pt")),
                #name = "cudcRNAseq",
                #column_split = factor(genepuree2$Growth, levels = c("Grow", "nogrow")),
                cluster_columns= FALSE,
                column_order=column_order(hm_z2),
                #row_order = row_order(hm_z2),
                show_column_dend = TRUE,
                cluster_column_slices = FALSE,
                
                cluster_rows=FALSE,
                
                #top_annotation=genepuree3,
                row_split = factor(gsva_u_sig$Set, levels = c(
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
                
                
                width = unit(2, "in"),  
                heatmap_height = unit(65, "in")
                
)



combo2<-hm_z2+hm_u2



gh2<-grid.grabExpr(draw(combo2))


ggsave("plots/hm-timex-all-rpkm-combo-clustered_relative-local.pdf",
       plot=gh2,
       
       scale = 1,
       dpi=600,
       width = 45,
       height = 100,
       unit="in",
       limitsize = FALSE
       
)


