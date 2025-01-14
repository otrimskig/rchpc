source("libs.R")
library(tidyverse)
library(dtplyr)
suppressPackageStartupMessages(library(ComplexHeatmap))
ht_opt$message = FALSE
library(gridtext)
library(circlize)
library(purrr)




pathways<-readRDS("acral_paired/ds/gsva_pathway_names.rds")%>%
  
  
  filter(list_name=="immune")%>%
  
  mutate(pathway_name_clean2=gsub("_", " ", pathway_name))




gsva_u<-readRDS("acral_paired/ds/gsva_u.rds")
gsva_z<-readRDS("acral_paired/ds/gsva_z.rds")

summary_stats<-readRDS("acral_paired/ds/gsva_pathway_stats.rds")

sample_info<-readRDS("acral_paired/ds/v00-sample_info.rds")





###now annotation set-up.

#this will be all the annotation information that we will want for 
#each plot. 

anno_df<-sample_info%>%select(sample_type)

anno_color_set <- list(
  sample_type = setNames(c("#a200ff", "#59f7c8"), 
                         c("acral", "subq"))
)



#get the names of the variables that we are interested in filtering. 
anno_colnames<-colnames(anno_df)


anno_subset<-list()
for (colnum in 1:length(anno_colnames)){
  
  #set column num(this wil be dynamic in for loop.)
  colnum<-1
  
  #get column name as character, from dynamically defined number.
  colname<-anno_colnames[colnum]
  
  #from that column name, get all the unique values.
  unique_to_filter<-anno_df%>%select(!!sym(colname))%>%unique()%>%pull()%>%as.character()


  
  anno_subset[[colname]]<-anno_color_set[[colname]][c(unique_to_filter)]
  
}


col_proper_names<-c("Tumor location")



anno<-HeatmapAnnotation(df=anno_df, 
                        col=c(anno_subset),
                        annotation_name_side = "left",
                        
                        gp = gpar(col = "black", fontsize = 3),
                        
                        annotation_label = c(col_proper_names),
                        annotation_name_gp= gpar(fontsize=7),
                        
                        simple_anno_size = unit(.125, "in"))



##########################################################

gene_lists<-unique(pathways$list_name)

# i<-1
#for (i in 1:length(gene_lists)){
# subset_factor<-gene_lists[i]

# pathways_to_include<-pathways%>%filter(list_name==subset_factor)%>%pull(gsva_name)%>%as.character()
pathways_to_include<-pathways%>%pull(gsva_name)%>%as.character()

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
  filter(pathway %in% c(pathways_to_include))%>%
  
  column_to_rownames("pathway")%>%
  as.matrix.data.frame()



row_labels<- structure(pathways$pathway_name_clean2, names = pathways$gsva_name)





mat_plus<-gsva_subset%>%
  as.data.frame()%>%
  rownames_to_column("pathway")%>%
  
  left_join(summary_stats_subset[,"mean_value"]%>%as.data.frame()%>%
              rownames_to_column("pathway")%>%
              rename("mean_value"=".")
  )%>%
  column_to_rownames("pathway")%>%
  as.matrix.data.frame()




###################################################################

ha3 = rowAnnotation(a=anno_points(mat[,"mean_value"],
                                  pch = c("|"), 
                                  #gp = gpar(col = 2:3))
                                  #joyplot_scale = 5,
                                  width = unit(2, "cm"),
                                  axis_param = list(
                                    side = "top")),
                    show_annotation_name = c(a = FALSE))
                                  

mat<-mat_plus[,-which(colnames(mat_plus) == "mean_value")]


hm1<- Heatmap(gsva_subset,
                right_annotation = ha3,
                 top_annotation = anno,
                column_title = gt_render(
                  paste0("<span style='font-size:25pt'>","Pathway Group: Immune ", "</span><br><span style='font-size:15pt'>expression levels normalized per row (pathway)</span>"), 
                  r = unit(2, "pt")),
              
                column_dend_height = unit(1, "in"),
                row_dend_width = unit(2, "in"),
                row_names_side = "right",
                row_gap = unit(1, "in"),
                row_labels = row_labels[rownames(gsva_subset)],
              
              
                
                column_title_gp = gpar(font = 2, fontsize = 60),
                row_title_gp = gpar(font = 2, fontsize = 40),
                show_heatmap_legend = FALSE,

                width = unit(10, "in"),
                heatmap_height = unit(nrow(gsva_subset)*.3, "in")
                
)



####################################



combo<-draw((hm1), annotation_legend_side = "left")

gh2<-grid.grabExpr(draw(combo))

ggsave(paste0("acral_paired/plots/hm-immune-clustered-relative-mean_exp_anno.pdf"),
       plot=gh2,
       
       
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/rchpc/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       
       
       scale = .8,
       dpi=600,
       width = 30,
       height = 25,
       unit="in",
       limitsize = FALSE
       
)




