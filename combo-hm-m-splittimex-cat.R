#combination heatmap with "z-score" heat. 

source("libs.R")
library(tidyverse)
library(dtplyr)
suppressPackageStartupMessages(library(ComplexHeatmap))
ht_opt$message = FALSE
library(gridtext)
library(circlize)
library(purrr)


#load relevant created from gsva analysis. 
# gsva_u_sig<-readRDS("ds/gsva_sig_u.rds")
gsva_z_sig<-readRDS("ds/gsva_sig_z.rds")
gsva_u<-readRDS("ds/gsva_u.rds")
gsva_z<-readRDS("ds/gsva_z.rds")


standard_error <- function(x) {
  sd(x) / sqrt(length(x))
}

summary_stats <- gsva_u%>%as.data.frame()%>%
  rownames_to_column("pathway")%>%
  rowwise() %>%
  mutate(
    min_value = min(c_across(`22925`:`36184`)),
    max_value = max(c_across(`22925`:`36184`)),
    mean_value = mean(c_across(`22925`:`36184`)),
    sd_value = sd(c_across(`22925`:`36184`)),
    se_value = standard_error(c_across(`22925`:`36184`))
  ) %>%
  select(-c(`22925`:`36184`))%>%
  ungroup()%>%
  column_to_rownames("pathway")%>%
  as.matrix.data.frame()

saveRDS(summary_stats, "ds/gsva_pathway_stats.rds")

anno_color<-readRDS("ds/colors_list.rds")
aod_colors = circlize::colorRamp2(c(50, 150), c("navy", "white"))


de_samples<-readRDS("ds/v10-per_sample_updated.rds")

###now annotation set-up.

#this will be all the annotation information that we will want for 
#each plot. note that this contains only the information relevant to 
#each specific plot.
de_anno_df<-de_samples%>%select(patho_grade, patho_cat_name, patho_cat2_name, patho_cat_det_name, aod, resultant_geno)%>%
  relocate(aod, .after = last_col())


#this is the df list that contains all assigned color values
#for all potential annotations.
#the issue is that it also contains colors for things 
#that WON'T show up in each plot.
#thus we need to subset it by only the variables we are interested in.
anno_color_set<-readRDS("ds/colors_list.rds")



#get the names of the variables that we are interested in filtering. 
de_anno_colnames<-colnames(de_anno_df)


anno_subset<-list()

for (colnum in 1:length(de_anno_colnames)){
  
  #set column num(this wil be dynamic in for loop.)
  # colnum<-1
  
  #get column name as character, from dynamically defined number.
  colname<-de_anno_colnames[colnum]
  
  #from that column name, get all the unique values.
  unique_to_filter<-de_anno_df%>%select(!!sym(colname))%>%unique()%>%pull()%>%as.character()
  
  
  
  
  anno_subset[[colname]]<-anno_color_set[[colname]][c(unique_to_filter)]
  
  
}



col_proper_names<-read_csv("ds/col_proper_names.csv")%>%filter(col_name %in% de_anno_colnames)%>%
  mutate(col_name=factor(col_name, levels=names(anno_subset)))%>%
  filter(!is.na(col_name))%>%arrange(col_name)%>%pull(proper_name)%>%as.character()




aod_colors = circlize::colorRamp2(c(50, 150), c("navy", "white"))


anno<-HeatmapAnnotation(df=de_anno_df, 
                        col=c(anno_subset, aod=aod_colors),
                        annotation_name_side = "left",
                        
                        gp = gpar(col = "black", fontsize = 3),
                        
                        annotation_label = c(col_proper_names, "Age of Death"),
                        annotation_name_gp= gpar(fontsize=7),
                        
                        simple_anno_size = unit(.125, "in"))



# anno<-HeatmapAnnotation(df=de_anno_df, 
#                         col=c(anno_color, aod=aod_colors),
#                         annotation_name_side = "left",
#                         gp = gpar(col = "black"),
#                         annotation_label = c(col_proper_names, "Age of Death"),
#                         annotation_name_gp= gpar(fontsize=7),
#                         simple_anno_size = unit(.25, "in"))

###########################################






















# 
# 
# 
# 
# hm_z <- Heatmap(gsva_z,
#                 
#                 column_dend_height = unit(1, "in"),
#                 row_dend_width = unit(2, "in"),
#                 
#                 
#                 top_annotation = anno,
#                 column_title = gt_render(
#                   paste0("<span style='font-size:25pt'>Expression Levels Normalized per gene <br>(Relative expression)</span>"), 
#                   r = unit(2, "pt")),
#                     #name = "cudcRNAseq",
#                     #column_split = factor(genepuree2$Growth, levels = c("Grow", "nogrow")),
#                     #cluster_columns= FALSE,
#                     cluster_column_slices = FALSE,
#                     #top_annotation=genepuree3,
#                     row_split = factor(gsva_z_sig$Set, levels = c(
#                       "Hallmark",
#                       "KEGG",
#                       "TIMEx",
#                       "Immune")),
#                 
#                 row_gap = unit(1, "in"),
#                 
#                 
#                 
#                 #sometimes fails not sure why 
#                     #cluster_row_slices = FALSE,
#                     #row_names_gp = gpar(fontsize = 6),
#                     #column_names_gp = gpar(fontsize =2),
#                     column_title_gp = gpar(font = 2, fontsize = 60),
#                     row_title_gp = gpar(font = 2, fontsize = 40),
#                     show_heatmap_legend = FALSE,
#                     
# 
#                  
#                     width = unit(10, "in"),
#                     heatmap_height = unit(65, "in"))
# 
# 
# hm_u <- Heatmap(gsva_u,
#                 column_title = gt_render(
#                   paste0("<span style='font-size:25pt'>Absolute <br>expression</span>"), 
#                   r = unit(2, "pt")),
#                 
#                 column_dend_height = unit(1, "in"),
#                 row_dend_width = unit(2, "in"),
#                 #col = colorRamp2(seq(min(gsva_u), max(gsva_u), length = 3), c("black", "white", "purple")),
#                 
#                 #name = "cudcRNAseq",
#                 #column_split = factor(genepuree2$Growth, levels = c("Grow", "nogrow")),
#                 cluster_columns= FALSE,
#                 column_order=column_order(hm_z),
#                 show_column_dend = TRUE,
#                 cluster_column_slices = FALSE,
#                 #top_annotation=genepuree3,
#                 row_split = factor(gsva_u_sig$Set, levels = c(
#                   "Hallmark",
#                   "KEGG",
#                   "TIMEx",
#                   "Immune")),
#                 
#                 
#                 
#                 
#                 
#                 row_gap = unit(1, "in"),
#                 
#                 
#                 
#                 
#                 #sometimes fails not sure why 
#                 #cluster_row_slices = FALSE,
#                 #row_names_gp = gpar(fontsize = 6),
#                 #column_names_gp = gpar(fontsize =2),
#                 column_title_gp = gpar(font = 2, fontsize = 60),
#                 row_title_gp = gpar(font = 2, fontsize = 40),
#                 show_heatmap_legend = FALSE,
#                 
#                 
#                 width = unit(1.5, "in"),  
#                 heatmap_height = unit(65, "in")
#                 
# )
# 
# str(hm_u)
# 
# combo<-draw((hm_u+hm_z),annotation_legend_side = "left")
# 
# 
# 
# gh<-grid.grabExpr(draw(combo))
# 
# 
# ggsave("plots/hm-timex-all-rpkm-combo-clustered_abs-local.pdf",
#        plot=gh,
#        
#        scale = 1,
#        dpi=600,
#        width = 45,
#        height = 100,
#        unit="in",
#        limitsize = FALSE
#        
# )


#######################################################################

gene_lists<-c("Hallmark",
                "KEGG",
                "TIMEx",
                "Immune")

# i<-1
for (i in 1:length(gene_lists)){
subset_factor<-gene_lists[i]

pathways_to_include<-gsva_z_sig%>%as.data.frame()%>%
  filter(Set==subset_factor)%>%pull(Signature)%>%as.character()


gsva_subset<-gsva_z%>%as.data.frame()%>%
  rownames_to_column("pathway")%>%
  filter(pathway %in% c(pathways_to_include))%>%
  column_to_rownames("pathway")%>%
  as.matrix.data.frame()

gsva_subset2<-gsva_u%>%as.data.frame()%>%
  rownames_to_column("pathway")%>%
  filter(pathway %in% c(pathways_to_include))%>%
  column_to_rownames("pathway")%>%
  as.matrix.data.frame()


summary_stats_subset<-summary_stats%>%as.data.frame()%>%
  rownames_to_column("pathway")%>%
  filter(pathway %in% c(pathways_to_include))%>%
  column_to_rownames("pathway")%>%
  as.matrix.data.frame()

summary_stats_subset2<-summary_stats_subset[,c(1,2,3)]


# ha = rowAnnotation(foo = anno_points(summary_stats_subset2,
#                                      width = unit(2, "cm"),
#                                      axis_param = list(
#                                        side = "bottom",
#                                        at = c(0, 0.5, 1), 
#                                        labels = c("zero", "half", "one"),
#                                        labels_rot = 45
#                                      ))
# )
# 
# ha2 = rowAnnotation(foo = anno_points(summary_stats_subset2[,2],
#                                       anno_simple(2, pch = 1), 
#                                       #gp = gpar(col = 2:3))
#                                        #joyplot_scale = 5,
#                                        width = unit(2, "cm"),
#                                      axis_param = list(
#                                        side = "bottom"
#                                   
#                                      ))
# )




# ha4 = rowAnnotation(foo = anno_density(gsva_subset2,
#                                        type="heatmap",
#                                        heatmap_colors = c("white","black"),
#                                       
#                                       width = unit(2, "cm"),
#                                       axis_param = list(
#                                         side = "bottom"
#                                       
#                                       )))

# 
# draw(ha3)
# draw(ha4)

ha3 = rowAnnotation(a=anno_points(summary_stats_subset[,3],
                                  pch = c("|"), 
                                  #gp = gpar(col = 2:3))
                                  #joyplot_scale = 5,
                                  width = unit(2, "cm"),
                                  axis_param = list(
                                    side = "top")),
                    show_annotation_name = c(a = FALSE))
                                  



hm1<- Heatmap(gsva_subset,
                right_annotation = ha3,
                 top_annotation = anno,
                column_title = gt_render(
                  paste0("<span style='font-size:25pt'>","Pathway Group: ", subset_factor, "</span><br><span style='font-size:15pt'>expression levels normalized per row (pathway)</span>"), 
                  r = unit(2, "pt")),
              
                column_dend_height = unit(1, "in"),
                row_dend_width = unit(2, "in"),
                row_names_side = "right",
                row_gap = unit(1, "in"),
                
                column_title_gp = gpar(font = 2, fontsize = 60),
                row_title_gp = gpar(font = 2, fontsize = 40),
                show_heatmap_legend = FALSE,

                width = unit(10, "in"),
                heatmap_height = unit(nrow(gsva_subset)*.3, "in")
                
)







########################################

# hm2 <- Heatmap(gsva_subset2,
#                 column_title = gt_render(
#                   paste0("<span style='font-size:25pt'>Absolute <br>expression</span>"),
#                   r = unit(2, "pt")),
#                
#                 cluster_columns= FALSE,
#                 column_order=column_order(hm_init),
#                
#                 show_column_dend = TRUE,
#                 cluster_column_slices = FALSE,
#                 cluster_rows=FALSE,
# 
#                 column_dend_height = unit(1, "in"),
#                 row_dend_width = unit(2, "in"),
# 
#                 row_gap = unit(1, "in"),
#                
#                 column_title_gp = gpar(font = 2, fontsize = 60),
#                 row_title_gp = gpar(font = 2, fontsize = 40),
#                 show_heatmap_legend = FALSE,
# 
#                 width = unit(1.5, "in"),
#                 heatmap_height = unit(nrow(gsva_subset)*.3, "in")
# )
# 


#####################################################



combo<-draw((hm1), annotation_legend_side = "left")

gh2<-grid.grabExpr(draw(combo))

ggsave(paste0("plots/hm-timex-", subset_factor, "-clustered-relative-summary-stats.pdf"),
       plot=gh2,
       
       scale = 1,
       dpi=600,
       width = 45,
       height = 50,
       unit="in",
       limitsize = FALSE
       
)



}
