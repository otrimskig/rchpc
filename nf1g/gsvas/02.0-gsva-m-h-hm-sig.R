source("libs.R")
suppressPackageStartupMessages(library(ComplexHeatmap))
ht_opt$message = FALSE
library(gridtext)
library(circlize)
library(purrr)
library(tidyverse)


gsva_stats<-readRDS("nf1g/gsvas/ds/gsva_pathway_stats_gliomas.rds")

gmatu0<-readRDS("nf1g/gsvas/ds/gsva_pathways_matrix.rds")

gsva_meta<-readRDS("nf1g/gsvas/ds/gsva_pathways_meta.rds")

sample_info<-readRDS("nf1g/ds/v10-per_sample_updated.rds")




n_row_total<-100

sig_pathways<-gsva_stats$results2%>%
  filter(!grepl("^m1", Pathway))%>%
  
  filter(min_pval<.05)%>%
  arrange(desc(abs_diff))%>%
  slice(1:n_row_total)
  




#subset sample metadata for selected samples.
samples_df<-sample_info%>%
  filter(patho_cat_name=="Gliomas")


#subset based on sample selection.
gmatu1<-gmatu0[sample(nrow(gmatu0), 200),samples_df$mouse_num]

gmatu1<-gmatu0[sig_pathways$Pathway ,samples_df$mouse_num]

#rescale matrix based on sample selection, rowwise.
gmatz1<-t(scale(t(gmatu1)))






anno_color_set<-readRDS("nf1g/ds/colors_list.rds")

de_anno_colnames<-c(#"patho_cat_name", 
                    "patho_grade", "resultant_geno")



anno_subset<-list()

for (colnum in 1:length(de_anno_colnames)){
  
  #set column num(this wil be dynamic in for loop.)
  # colnum<-1
  
  #get column name as character, from dynamically defined number.
  colname<-de_anno_colnames[colnum]
  
  #from that column name, get all the unique values.
  unique_to_filter<-samples_df%>%select(!!sym(colname))%>%unique()%>%pull()%>%as.character()
  
  anno_subset[[colname]]<-anno_color_set[[colname]][c(unique_to_filter)]
  
}






anno<-HeatmapAnnotation(df=samples_df%>%
                          select(mouse_num, 
                                 patho_grade,
                                 resultant_geno), 
                        col=c(anno_subset),
                        annotation_name_side = "left",
                        
                        gp = gpar(col = "black", fontsize = 3),
                        
                        #annotation_label = c(col_proper_names),
                        annotation_name_gp= gpar(fontsize=7),
                        
                        simple_anno_size = unit(.125, "in"))



size <- dev.size("in")
scaled_height <- unit(size[2] * 2, "in")  # 90% of available height
scaled_width <- unit(size[1] * 0.9, "in")   # 90% of available width


hm_print <- Heatmap(gmatz1,
                    top_annotation = anno,
                    column_dend_height = unit(1, "in"),
                    row_dend_width = unit(2, "in"),
                    row_names_side = "right",
                    row_gap = unit(1, "in"),
                    column_title_gp = gpar(font = 2, fontsize = 60),
                    row_title_gp = gpar(font = 2, fontsize = 1),
                    
                    row_names_gp = gpar(fontsize = 8),
                    
                    show_heatmap_legend = FALSE,
                    width = scaled_width,
                    heatmap_height = scaled_height)



draw(hm_print, 
     heatmap_legend_side = "top", 
     annotation_legend_side = "left")














combo<-draw((hm_u+hm_z),annotation_legend_side = "left")

combo<-hm1

gh<-grid.grabExpr(draw(combo))


ggsave("nf1g/gsvas/plots/hm-gsva-500.pdf",
       plot=gh,
       scale = 1,
       dpi=600,
       width = 45,
       height = 100,
       unit="in",
       limitsize = FALSE
)

