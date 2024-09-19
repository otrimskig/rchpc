source("libs.R")
library(tidyverse)
library(dtplyr)
library(gridtext)
library(purrr)

suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
ht_opt$message = FALSE

#custom functions split_matrix and standard error.
split_matrix <- function(matrix, segment_size) {
  
  # Create a sequence of indices for each segment
  split_indices <- split(seq_len(num_rows), ceiling(seq_len(num_rows) / segment_size))
  
  # Use these indices to split the matrix into a list of data frames
  split_list <- lapply(split_indices, function(indices) matrix[indices, , drop = FALSE])
  
  return(split_list)
}


standard_error <- function(x) {
  sd(x) / sqrt(length(x))
}


#datasets to be used in all loops.
gsva_analysed<-readRDS("nf1g/ds/gsva_analysis_ds_list_ms.rds")
de_samples<-readRDS("nf1g/ds/v10-per_sample_updated.rds")

####################################
#cut to check loop. 
gsva_analysed<-gsva_analysed[1]
  
  
for (i in 1:length(gsva_analysed)){
  
  
#PART 1: get in unscaled matrix and clean. Calculate scaled matrix. 
  
#get name of analysis performed (which large signature list used.)  
analysis_name<-names(gsva_analysed[i])


#get unscaled signature scores for all samples.  
gsva_u<-readRDS(gsva_analysed[[i]][["u"]][[1]])
  
#some quick fixing to ensure sample_ids are changed to mouse_nums
#and to only include proper samples. 
old_names<-de_samples$sample_id
new_names<-de_samples$mouse_num
name_mapping <- setNames(new_names, old_names)
  
# check name lengths.
length(old_names)==length(new_names)
  
# filter matrix for nf1g samples that should be kept,
#and change from sampleids to mouse_nums, for gsva_u. 
gsva_u<-gsva_u[, colnames(gsva_u) %in% old_names]
colnames(gsva_u) <- name_mapping[colnames(gsva_u)]
  
  
#after cleaning unsclaed matrix and excluding proper samples, THEN scale by row. 
gsva_z<-t(scale(t(gsva_u)))


  
#PART 2a: calculate stats for scaled (z-score) df. Exclude ones with smaller differences. (filter(diff>=..))

#get name of columns to stat on
columns_to_select_for_stats1<-colnames(gsva_z)

#find difference between max and min values,
#subset based on calculated differences.
large_diffs<-gsva_z%>%as.data.frame()%>%
    rownames_to_column("pathway")%>%
    rowwise()%>%
    mutate(min_value = min(across(all_of(columns_to_select_for_stats1))),
           max_value = max(across(all_of(columns_to_select_for_stats1)))
    )%>%
    ungroup()%>%
    mutate(diff=max_value-min_value)%>%
    
  
  #filter(diff>=5)%>%
    
  pull("pathway")
  
#sort rownames in calculated dataset alphabetically,
#and subset by calculated threshold.
uro<-gsva_u[order(rownames(gsva_u)), ][large_diffs,]
zro<-gsva_z[order(rownames(gsva_z)), ][large_diffs,]



#PART 2b: calculate stats for unscaled df. A subset of this will be used to make the side plot on the hms.

#make a summary dataframe 
summary_stats_u<-gsva_u%>%as.data.frame()%>%
  rownames_to_column("pathway")%>%
  rowwise()

#you have to redefine the column subset - otherwise mean() won't work. 
columns_to_select_for_stats2 <- colnames(summary_stats_u)[colnames(summary_stats_u) != "pathway"]
  
summary_stats_u<-summary_stats_u%>%
  mutate(min_value = min(across(all_of(columns_to_select_for_stats2))),
         max_value = max(across(all_of(columns_to_select_for_stats2))))%>%
  
  rowwise() %>%
  mutate(
    mean_value = mean(c_across(all_of(columns_to_select_for_stats2)), na.rm = TRUE),
    sd_value = sd(c_across(all_of(columns_to_select_for_stats2)), na.rm = TRUE),
         sd_value = sd(across(all_of(columns_to_select_for_stats2))),
         se_value = standard_error(across(all_of(columns_to_select_for_stats2)))
  )%>%
  select(-all_of(columns_to_select_for_stats1))%>%
  ungroup()%>%
  column_to_rownames("pathway")%>%
  as.matrix.data.frame()






#PART 3: Split matrix into manageable chunks. 
#Based on that, subset all necessary data and pass to build hms.

#first check that the dataset after filtering still retains data.
#otherwise, building heatmap will fail. 



if (is.null(nrow(uro))==FALSE){

num_rows<-nrow(uro)
#creates splitting method.
split_segments <- split_matrix(uro, 1999)
#gets the numbers (how many split exist). Use for indexing in passing to loop.
segments_nums<-names(split_segments)


###################################
#remove when doing loop.s
# segments_nums<-segments_nums[1]




for (seg in 1:length(segments_nums)){

seg_num<-sprintf("%03d", seg)  
  
  
#based on split matrix, get signatures to keep.
keep_signatures<-rownames(split_segments[[seg]])
  
#subset the data based on keeping signatures. 
gsva_u_seg<-uro[keep_signatures,]
gsva_z_seg<-zro[keep_signatures,]

pathways_to_include<-rownames(gsva_u_seg)

gsva_subset<-gsva_z%>%as.data.frame()%>%
  rownames_to_column("pathway")%>%
  filter(pathway %in% c(pathways_to_include))%>%
  column_to_rownames("pathway")%>%
  as.matrix.data.frame()

summary_stats_subset<-summary_stats_u%>%as.data.frame()%>%
  rownames_to_column("pathway")%>%
  filter(pathway %in% c(pathways_to_include))%>%
  column_to_rownames("pathway")%>%
  as.matrix.data.frame()



###now annotation set-up.

#this will be all the annotation information that we will want for 
#each plot. note that this contains only the information relevant to 
#each specific plot.
de_anno_df<-de_samples%>%select(patho_grade, patho_cat_name, 
                                # patho_cat_det_name, 
                                resultant_geno)


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





ha3 = rowAnnotation(a=anno_points(summary_stats_subset[,3],
                                  pch = c("|"), 
                                  #gp = gpar(col = 2:3))
                                  #joyplot_scale = 5,
                                  width = unit(2, "cm"),
                                  axis_param = list(
                                    side = "top")),
                    show_annotation_name = c(a = FALSE))
                                  

names2<-de_samples%>%
  select(patho_cat_name)%>%
  group_by(patho_cat_name)%>%slice(1)%>%ungroup()%>%
  pull()


nums2<-de_samples%>%
  filter(patho_cat_name==names2[1])%>%
  pull(mouse_num)


mat1<-gsva_subset[,nums2]
  

hm1<- Heatmap(mat1,
              #right_annotation = ha3,
              #top_annotation = anno,
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
              #width = unit(10, "in"),
              heatmap_height = unit(nrow(gsva_subset)*.3+5, "in")
              
             
              
)





nums2<-de_samples%>%
  filter(patho_cat_name==names2[2])%>%
  pull(mouse_num)


mat2<-gsva_subset[,nums2]



hm2<- Heatmap(mat2,
              #right_annotation = ha3,
              #top_annotation = anno,
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
              #width = unit(10, "in"),
              heatmap_height = unit(nrow(gsva_subset)*.3+5, "in")
              
)




nums2<-de_samples%>%
  filter(patho_cat_name==names2[3])%>%
  pull(mouse_num)


mat3<-gsva_subset[,nums2]



hm3<- Heatmap(mat3,
              #right_annotation = ha3,
              #top_annotation = anno,
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
              #width = unit(10, "in"),
              heatmap_height = unit(nrow(gsva_subset)*.3+5, "in")
              
            
              
)





nums2<-de_samples%>%
  filter(patho_cat_name==names2[4])%>%
  pull(mouse_num)


mat4<-gsva_subset[,nums2]



hm4<- Heatmap(mat4,
              #right_annotation = ha3,
              #top_annotation = anno,
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
              #width = unit(10, "in"),
              heatmap_height = unit(nrow(gsva_subset)*.3+5, "in")
              
              
              
)




nums2<-de_samples%>%
  filter(patho_cat_name==names2[5])%>%
  pull(mouse_num)


mat5<-gsva_subset[,nums2]



hm5<- Heatmap(mat5,
              right_annotation = ha3,
              #top_annotation = anno,
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
              #width = unit(10, "in"),
              heatmap_height = unit(nrow(gsva_subset)*.3+5, "in")
              
              
              
)








hm_split<-hm1+hm2+hm3+hm4+hm5



#####################################################


hm_name_pdf<-paste0("nf1g/plots/gsva_hms/hm-largeplots-disease-groups", analysis_name, "-segment-", seg_num, ".pdf")


gh2<-grid.grabExpr(draw(draw(hm_split, 
                             ht_gap = unit(10, "mm"),
                             heatmap_width=unit(10, "in"),
                             annotation_legend_side = "left")))

ggsave(hm_name_pdf,
       
       plot=gh2,
       
       scale = 1,
       dpi=600,
       width = 20,
       height = 500,
       unit="in",
       limitsize = FALSE
       
)

#status update
print(paste0(round(Sys.time()), " : ", hm_name_pdf, " written"))



# combo<-draw((hm1), annotation_legend_side = "left")

# gh2<-grid.grabExpr(draw(draw((hm_split), annotation_legend_side = "left")))
# 
# ggsave(sub(".pdf$", "-rows_unclustered.pdf", hm_name_pdf), 
#        
#        plot=gh2,
#        
#        scale = 1,
#        dpi=600,
#        width = 45,
#        height = 600,
#        unit="in",
#        limitsize = FALSE
#        
# )

#status update
print(paste0(round(Sys.time()), " : ", sub(".pdf$", "-rows_unclustered.pdf", hm_name_pdf), " written"))

  
  
}

}else{
  
print("no obs.")  

  }
}
  
  
  
stop("end of maker loops. to compile pdfs finish code.")



plots<-fs::dir_info("nf1g/plots/gsva_hms", full.names = T, pattern="^hm.*\\.pdf$", recurse = F)%>%tibble()%>%
  filter(modification_time>Sys.time()-lubridate::hours(4))%>%
  pull(path)



qpdf::pdf_combine(plots, output = "nf1g/plots/combined/hm-gsva-combined-2.pdf")
