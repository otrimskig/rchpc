#combination heatmap with "z-score" heat. 
source("libs.R")
library(dtplyr)
suppressPackageStartupMessages(library(ComplexHeatmap))
ht_opt$message = FALSE
library(gridtext)
library(circlize)
library(purrr)
standard_error <- function(x) {
  sd(x) / sqrt(length(x))
}
library(tidyverse)


#load relevant created from gsva analysis. 
sample_info<-readRDS("nf1g/ds/v10-per_sample_updated.rds")

gsva_pathway_stats<-readRDS("nf1g/gsvas/gsva_top500_output.rds")$results2%>%
  arrange(desc(abs_diff))%>%
  slice(1:100)

pathway_meta<-readRDS("nf1g/ds/gsva/gsva_pathways_meta.rds")



gsva_u0<-readRDS("nf1g/ds/gsva/gsva_pathways_matrix.rds")[gsva_pathway_stats$Pathway,]

samples<-sample_info%>%
  filter(patho_cat_name=="Gliomas")%>%
  pull(sample_id)

gsva_u1<-gsva_u0[,samples]

gsva_z1<-t(scale(t(gsva_u1)))


###now annotation set-up.
#this will be all the annotation information that we will want for 
#each plot. note that this contains only the information relevant to 
#each specific plot.
de_anno_df<-sample_info%>%
  filter(patho_cat_name=="Gliomas")%>%
  dplyr::select(patho_grade, patho_cat_det_name, resultant_geno)



#this is the df list that contains all assigned color values
#for all potential annotations.
#the issue is that it also contains colors for things 
#that WON'T show up in each plot.
#thus we need to subset it by only the variables we are interested in.
anno_color_set<-readRDS("nf1g/ds/colors_list.rds")



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


anno<-HeatmapAnnotation(df=de_anno_df, 
                        col=c(anno_subset),
                        annotation_name_side = "left",
                        
                        gp = gpar(col = "black", fontsize = 3),
                        
                        #annotation_label = c(col_proper_names),
                        annotation_name_gp= gpar(fontsize=7),
                        
                        simple_anno_size = unit(.125, "in"))


hm1<- Heatmap(gsva_z1,
              #right_annotation = ha3,
              top_annotation = anno,
              #column_title = gt_render(
              # paste0("<span style='font-size:25pt'>","Pathway Group: ", subset_factor, "</span><br><span style='font-size:15pt'>expression levels normalized per row (pathway)</span>"), 
              #  r = unit(2, "pt")),
              
              
              column_dend_height = unit(1, "in"),
              row_dend_width = unit(2, "in"),
              row_names_side = "right",
              row_gap = unit(1, "in"),
              
              column_title_gp = gpar(font = 2, fontsize = 60),
              row_title_gp = gpar(font = 2, fontsize = 40),
              show_heatmap_legend = FALSE,
              
              width = unit(3, "in"),
              heatmap_height = unit(nrow(gsva_z1)*.3, "in"))

hm1

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


