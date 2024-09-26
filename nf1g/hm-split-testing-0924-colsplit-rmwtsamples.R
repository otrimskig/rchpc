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
# gsva_analysed<-gsva_analysed[1]
# 
# i<-1

for (i in 1:length(gsva_analysed)){



#PART 1: get in unscaled matrix and clean. Calculate scaled matrix. 

# #get name of analysis performed (which large signature list used.)  
# analysis_name<-names(gsva_analysed[i])



#get unscaled signature scores for all samples.  
g0u<-readRDS(gsva_analysed[[i]][["u"]][[1]])

# #sort rownames in calculated dataset alphabetically
g1u<-g0u[order(rownames(g0u)), ]


#some quick fixing to ensure sample_ids are changed to mouse_nums
#get sample_ids and corresponding mouse_nums, for those included in
#each specific de_samples dataset.
#maps sample_ids and mouse_nums together. 
#sample_ids are the names, and mouse_nums are the values. 
sample_filter <- setNames(de_samples$mouse_num, de_samples$sample_id)

# filter matrix for samples that should be kept,
g2u<-g1u[, colnames(g1u) %in% names(sample_filter)]

#now change from sample ids to mouse_nums. 
colnames(g2u)<-sample_filter[colnames(g2u)]



#AFTER cleaning unscaled matrix and excluding proper samples, THEN scale by row. 
g2z<-t(scale(t(g2u)))
  
  
  





#PART 2-a: calculate stats for scaled (z-score) df. Exclude ones with smaller differences. (filter(diff>=..))

#get name of columns to stat on
columns_to_select_for_stats1<-colnames(g2z)


#find difference between max and min values,
#subset based on calculated differences.
diff_stats<-g2z%>%as.data.frame()%>%
  rownames_to_column("pathway")%>%
  rowwise()%>%
  mutate(min_value = min(across(all_of(columns_to_select_for_stats1))),
         max_value = max(across(all_of(columns_to_select_for_stats1)))
  )%>%
  ungroup()%>%
  mutate(diff=max_value-min_value)



#PART 2-b: calculate stats for unscaled df. A subset of this will be used to make the side plot on the hms.

#make a summary dataframe 
summary_stats_u<-g2u%>%as.data.frame()%>%
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



#part 2-c : subset pathways based on whatever criteria.
pathway_selection<-diff_stats%>%
  
  # filter(diff>=4)%>%
  
  pull("pathway")


#filter dataset.
g3u<-g2u[pathway_selection,]
g3z<-g2z[pathway_selection,]






#PART 3: Split matrix into manageable chunks. 
#Based on that, subset all necessary data and pass to build hms.

#first check that the dataset after filtering still retains data.
#otherwise, building heatmap will fail. 



if (is.null(nrow(g3u))==FALSE){
  
num_rows<-nrow(g3z)
  #creates splitting method.
split_segments <- split_matrix(g3z, 1999)
  #gets the numbers (how many split exist). Use for indexing in passing to loop.
segments_nums<-names(split_segments)
  
  
###################################
#remove when doing loop.s
# segments_nums<-segments_nums[1]
  
for (seg in 1:length(segments_nums)){
  
# seg<-1

seg_num<-sprintf("%03d", seg)  

#based on split matrix, get signatures to keep.
keep_signatures<-rownames(split_segments[[seg]])

#subset the data based on keeping signatures. 
gu_sub1<-g3u[keep_signatures,]
gz_sub1<-g3z[keep_signatures,]
summary_stats_u_sub1<-summary_stats_u[keep_signatures,]
pathways_to_include<-rownames(gu_sub1)
sample_names<-colnames(gu_sub1)




#this will be all the annotation information that we will want for 
#each plot. note that this contains only the information relevant to 
#each specific plot.
anno_df1<-de_samples%>%
  
  filter(sample_names %in% mouse_num)%>%
  relocate(patho_cat_name, patho_grade, resultant_geno)%>%
  
  select(mouse_num, 
         patho_cat_name, 
         patho_grade,
        # patho_cat_det_name, 
        resultant_geno)
  
  


#this is the df list that contains all assigned color values
#for all potential annotations.
#the issue is that it also contains colors for things 
#that WON'T show up in each plot.
#thus we need to subset it by only the variables we are interested in.
anno_color_set<-readRDS("nf1g/ds/colors_list.rds")



#get the names of the variables that we are interested in filtering. 
anno_info<-anno_df1%>%select(-mouse_num)


anno_subset<-list()

# 
# colnames(anno_info)


for (colnum in 1:length(colnames(anno_info))){
  
  #get column name as character, from dynamically defined number.
  colname<-colnames(anno_info)[colnum]
  
  #from that column name, get all the unique values.
  unique_to_filter<-anno_info%>%select(!!sym(colname))%>%unique()%>%pull()%>%as.character()
  
  anno_subset[[colname]]<-anno_color_set[[colname]][c(unique_to_filter)]
  
}


col_proper_names<-suppressMessages(read_csv("nf1g/ds/col_proper_names.csv"))%>%filter(col_name %in% colnames(anno_info))%>%
  mutate(col_name=factor(col_name, levels=names(anno_subset)))%>%
  filter(!is.na(col_name))%>%arrange(col_name)%>%pull(proper_name)%>%as.character()







column_subset<-anno_df1%>%group_by(patho_cat_name)%>%
  slice(1)%>%
  ungroup()%>%
  pull(patho_cat_name)



anno_block<-anno_df1

gz_sub1_block<-gz_sub1[,pull(anno_block, mouse_num)]


  
anno<-HeatmapAnnotation(df=anno_block%>%select(-mouse_num),
                         
                          col=c(anno_subset),
                          annotation_name_side = "left",
                          
                          
                          gp = gpar(col = "black", fontsize = 3),
                          
                          annotation_label = c(col_proper_names),
                          annotation_name_gp= gpar(font=1,fontsize=8),
                          
                          simple_anno_size = unit(4, "mm"))
  
 
ha1 = rowAnnotation(a=anno_points(summary_stats_u_sub1[,3],
                                  pch = c("|"), 
                                  gp = gpar(font = 4, fontsize = 8),
                                  #joyplot_scale = 5,
                                  
                                  
                                  
                                  width = unit(10, "mm"),
                                  
                                  axis_param = list(
                                    side = "top")),
                    show_annotation_name = c(a = FALSE))



assign(paste0("hm_1"),Heatmap(gz_sub1_block,
              #right_annotation = ha3,
              top_annotation = anno,
              
              column_split = anno_block$patho_cat_name,
              # column_title_rot = 90,
              column_title = NULL,
              
              column_gap = unit(4, "mm"),
              
              # column_title = gt_render(
              #   paste0("<span style='font-size:25pt'>","Pathway Group: ", analysis_name, " subset", seg_num, "</span><br><span style='font-size:15pt'>expression levels normalized per row (pathway)</span>"), 
              #   r = unit(2, "pt")),
              # 
              column_dend_height = unit(12, "mm"),
              row_dend_width = unit(16, "mm"),
              row_names_side = "right",
              
              row_names_gp = gpar(font = 1, fontsize = 8),
              
              
              row_gap = unit(1, "mm"),
              
              
              
              
              
             right_annotation =ha1,
              
              height = unit(nrow(gz_sub1_block)*3, "mm"),
              width = unit(ncol(gz_sub1_block)*8, "mm"),
              
              column_title_gp = gpar(font = 2, fontsize = 60),
              
              show_heatmap_legend = FALSE
              
              
             
              
))



# 






comb_hm<-hm_1
  
  # rowAnnotation(rn = anno_text(rownames(gz_sub1_block), 
  #                                                                              gp=gpar(font = 1, fontsize = 8),
  #                                                                              location = unit(0, "npc"), just = "left"))

# draw(comb_hm, annotation_legend_side = "left")











analysis_name<-names(gsva_analysed[i])

analysis_name
seg_num

hm_pdf_name<-paste0("nf1g/plots/gsva_hms/hm-disease-groups-clus-", analysis_name, "-segment-", seg_num, ".pdf")



hm_object_export<-grid.grabExpr(draw(draw(comb_hm, annotation_legend_side = "left")))

ggsave(hm_pdf_name,
       
       plot=hm_object_export,
       
       scale = 1,
       dpi=600,
       width = 700,
       height=nrow(gz_sub1_block)*3+200,
       
      
       unit="mm",
       limitsize = FALSE
       
)


print(paste0(round(Sys.time()), " : ", hm_pdf_name, " written"))
 
  
}
}else{
  
  print("no obs.") 
}
}





plots<-fs::dir_info("nf1g/plots/gsva_hms", full.names = T, pattern="^hm.*\\.pdf$", recurse = F)%>%tibble()%>%
  filter(modification_time>Sys.time()-lubridate::hours(4))%>%
  pull(path)



qpdf::pdf_combine(plots, output = "nf1g/plots/combined/hm-disease-group01.pdf")








