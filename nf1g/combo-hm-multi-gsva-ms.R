source("libs.R")
library(tidyverse)
library(dtplyr)
library(gridtext)
library(purrr)

suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
ht_opt$message = FALSE
# 
# if (!exists("n.cores")) {
#   # Run your code here
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
# }


# foreach(i=1:length(gsva_analysed)) %dopar% {





gsva_analysed<-readRDS("nf1g/ds/gsva_analysis_ds_list_ms.rds")

de_samples<-readRDS("nf1g/ds/v10-per_sample_updated.rds")


  
 
  
  
  
  
for (i in 1:length(gsva_analysed)){
  
  
  
  
analysis_name<-names(gsva_analysed[i])

  
  gsva_u<-readRDS(gsva_analysed[[i]][["u"]][[1]])
  gsva_z<-readRDS(gsva_analysed[[i]][["z"]][[1]])
  
  
  uro<-gsva_u[order(rownames(gsva_u)), ]
  zro<-gsva_z[order(rownames(gsva_z)), ]
  
  num_rows<-nrow(uro)

  
  
if (num_rows>0){
  


  
  split_matrix <- function(matrix, segment_size) {
    # Create a sequence of indices for each segment
    split_indices <- split(seq_len(num_rows), ceiling(seq_len(num_rows) / segment_size))
    
    # Use these indices to split the matrix into a list of data frames
    split_list <- lapply(split_indices, function(indices) matrix[indices, , drop = FALSE])
    
    return(split_list)
  }




split_segments <- split_matrix(uro, 100)
segments_nums<-names(split_segments)


for (seg in 1:length(segments_nums)){

seg_num<-sprintf("%03d", seg)
keep_signatures<-rownames(split_segments[[seg]])
  
gsva_u<-uro[keep_signatures,]
gsva_z<-zro[keep_signatures,]

hm_name_pdf<-paste0("nf1g/plots/gsva_hms/hm-", analysis_name, "segment-", seg_num, ".pdf")

old_names<-de_samples$sample_id
new_names<-de_samples$mouse_num
name_mapping <- setNames(new_names, old_names)

# check name lengths.
length(old_names)==length(new_names)


# filter matrix for nf1g samples that should be kept,
#and change from sampleids to mouse_nums, for 
#gsva_u. 
gsva_u<-gsva_u[, colnames(gsva_u) %in% old_names]
colnames(gsva_u) <- name_mapping[colnames(gsva_u)]

#repeat for gsva_z
gsva_z<-gsva_z[, colnames(gsva_z) %in% old_names]
colnames(gsva_z) <- name_mapping[colnames(gsva_z)]





standard_error <- function(x) {
  sd(x) / sqrt(length(x))
}

summary_stats_a <- gsva_u%>%as.data.frame()%>%
  rownames_to_column("pathway")%>%
  rowwise()

columns_to_select_for_stats<-colnames(summary_stats_a)[2:length(colnames(summary_stats_a))]

summary_stats<-summary_stats_a%>%
  mutate(min_value = min(c_across(all_of(columns_to_select_for_stats))),
         max_value = max(c_across(all_of(columns_to_select_for_stats))),
         mean_value = mean(c_across(all_of(columns_to_select_for_stats))),
         sd_value = sd(c_across(all_of(columns_to_select_for_stats))),
         se_value = standard_error(c_across(all_of(columns_to_select_for_stats)))
  )%>%
  select(-all_of(columns_to_select_for_stats))%>%
  ungroup()%>%
  column_to_rownames("pathway")%>%
  as.matrix.data.frame()







###now annotation set-up.

#this will be all the annotation information that we will want for 
#each plot. note that this contains only the information relevant to 
#each specific plot.
de_anno_df<-de_samples%>%select(patho_grade, patho_cat_name, 
                                #patho_cat2_name,
                                patho_cat_det_name, resultant_geno)


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
  
  #get column name as character, from dynamically defined number.
  colname<-de_anno_colnames[colnum]
  
  #from that column name, get all the unique values.
  unique_to_filter<-de_anno_df%>%select(!!sym(colname))%>%unique()%>%pull()%>%as.character()

  anno_subset[[colname]]<-anno_color_set[[colname]][c(unique_to_filter)]
  
}



col_proper_names<-suppressMessages(read_csv("nf1g/ds/col_proper_names.csv"))%>%filter(col_name %in% de_anno_colnames)%>%
  mutate(col_name=factor(col_name, levels=names(anno_subset)))%>%
  filter(!is.na(col_name))%>%arrange(col_name)%>%pull(proper_name)%>%as.character()



anno<-HeatmapAnnotation(df=de_anno_df, 
                        col=c(anno_subset),
                        annotation_name_side = "left",
                        
                        gp = gpar(col = "black", fontsize = 3),
                        
                        annotation_label = c(col_proper_names),
                        annotation_name_gp= gpar(fontsize=7),
                        
                        simple_anno_size = unit(.125, "in"))

#######################################################################


pathways_to_include<-rownames(gsva_u)

gsva_subset<-gsva_z%>%as.data.frame()%>%
  rownames_to_column("pathway")%>%
  filter(pathway %in% c(pathways_to_include))%>%
  column_to_rownames("pathway")%>%
  as.matrix.data.frame()


summary_stats_subset<-summary_stats%>%as.data.frame()%>%
  rownames_to_column("pathway")%>%
  filter(pathway %in% c(pathways_to_include))%>%
  column_to_rownames("pathway")%>%
  as.matrix.data.frame()



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
                  paste0("<span style='font-size:25pt'>","Pathway Group: ", analysis_name, " subset", seg_num, "</span><br><span style='font-size:15pt'>expression levels normalized per row (pathway)</span>"), 
                  r = unit(2, "pt")),
              
                column_dend_height = unit(1, "in"),
                row_dend_width = unit(2, "in"),
                row_names_side = "right",
                row_gap = unit(1, "in"),
                
                column_title_gp = gpar(font = 2, fontsize = 60),
                row_title_gp = gpar(font = 2, fontsize = 40),
                show_heatmap_legend = FALSE,

                width = unit(10, "in"),
                heatmap_height = unit(nrow(gsva_subset)*.3+1, "in")
                
)




hm2<- Heatmap(gsva_subset,
              right_annotation = ha3,
              top_annotation = anno,
              column_title = gt_render(
                paste0("<span style='font-size:25pt'>","Pathway Group: ", analysis_name, " subset", seg_num, "</span><br><span style='font-size:15pt'>expression levels normalized per row (pathway)</span>"), 
                r = unit(2, "pt")),
              
              column_dend_height = unit(1, "in"),
              
              row_order = sort(rownames(gsva_u)), 
              
              
              
              row_names_side = "right",
              row_gap = unit(1, "in"),
              
              column_title_gp = gpar(font = 2, fontsize = 60),
              row_title_gp = gpar(font = 2, fontsize = 40),
              show_heatmap_legend = FALSE,
              
              width = unit(10, "in"),
              heatmap_height = unit(nrow(gsva_subset)*.3+1, "in")
              
)











#####################################################

# combo<-draw((hm1), annotation_legend_side = "left")

gh2<-grid.grabExpr(draw(draw((hm1), annotation_legend_side = "left")))

ggsave(hm_name_pdf,
       
       plot=gh2,
       
       scale = 1,
       dpi=600,
       width = 45,
       height = 50,
       unit="in",
       limitsize = FALSE
       
)




# combo<-draw((hm1), annotation_legend_side = "left")

gh2<-grid.grabExpr(draw(draw((hm2), annotation_legend_side = "left")))

ggsave(sub(".pdf$", "-rows_unclustered.pdf", hm_name_pdf), 
       
       plot=gh2,
       
       scale = 1,
       dpi=600,
       width = 45,
       height = 50,
       unit="in",
       limitsize = FALSE
       
)




  
  
}
}else{
  
print("no obs.")  
}
}
  
  
  

library(tidyverse)


plots<-fs::dir_info("nf1g/plots/gsva_hms", full.names = T, pattern="^hm.*\\.pdf$", recurse = F)%>%tibble()%>%
  filter(modification_time>Sys.time()-lubridate::hours(4))%>%
  pull(path)



qpdf::pdf_combine(plots, output = "nf1g/plots/combined/hm-gsva-combined-1.pdf")