if (!exists("n.cores")) {
  # Run your code here
"initilizing cores..."
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

"parallel cores initialized."
}

source("libs.R")
library(tidyverse)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(circlize)
library(foreach)

###########################################
#set threshold values for deexps datasets.

logfc_threshold<-2
fdr_threshold<-.001

#get vector of all .rds files containing dexp analyses.
all_dexps<-list.files("m-dexps", full.names = T, pattern="^dexp.*\\.rds$", recursive = F)



#get sample metadata.
all_sample_info<-readRDS("ds/v10-per_sample_updated.rds")

##############################################

#start parallelized loop to make plot for each deexp dataset.

foreach(i=1:length(all_dexps)) %dopar% {

#for parallel loops, each library needs to be initialized inside the loop.  
source("libs.R") 
 library(tidyverse)
  library(ComplexHeatmap)
#   
# i<-3

rds_file_path<-all_dexps[i]






#determine plot title from basic naming pattern of rds file.
plot_title<-basename(rds_file_path)%>%
  sub(".rds$", "", .)%>%
  sub("^dexp-", "", .)%>%
  sub(".*-", "", .)%>%
  sub(" v. ", " vs. ", .)


#get file base since we will need this later to actually name the file. 
file_base<-basename(rds_file_path)%>%
  sub(".rds$", "", .)%>%
  sub("^dexp-", "", .)



#now load the actual dataset from the predetermined path.
de_df<-readRDS(rds_file_path)

#files have columns with rpkm_[sample_id] and the corresponding rpkm values.
#we can use that when we make the matrix for the heatmap - 
#but ultimately we want the mouse number as the column labels.
#get a vector of the (loaded) column names.
sa<-colnames(de_df)

#now use sample metadata to match those column names and create
#a vector of mouse_nums in the same order. Will be used for 
#labels in the columns of the heatmap.

column_labels_mouse_num<-tibble(col_name=sa[grep("^rpkm", sa)])%>%
  mutate(sample_id=sub("rpkm_", "", col_name))%>%
  left_join(all_sample_info%>%select(sample_id, mouse_num))%>%
  arrange(col_name)%>%
  pull(mouse_num)


#get the remaining metadata for those samples
de_samples<-tibble(mouse_num=column_labels_mouse_num)%>%
  left_join(all_sample_info)%>%
  mutate(patho_grade=as_factor(patho_grade))
  # mutate(patho_cat2=as_factor(patho_cat2))






hm_mat<-de_df%>%
  filter(abs(logFC)>=logfc_threshold)%>%
  filter(FDR<fdr_threshold)%>%
  
  select(-c(gene_id_ms:FDR))%>%
  column_to_rownames("gene_name_ms")



if (nrow(hm_mat) > 500) {
  # If true, sample 500 rows
  hm_mat_500 <- hm_mat %>% sample_n(500)
} else {
  # If false, don't sample
  hm_mat_500 <- hm_mat
}


scaled_mat<-t(scale(t(hm_mat_500)))



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
                        annotation_name_side = "right",
                        gp = gpar(col = "black", fontsize = 3),
                       
                        annotation_label = c(col_proper_names, "Age of Death"),
                        annotation_name_gp= gpar(fontsize=7),
                        
                        simple_anno_size = unit(.125, "in"))



width_scale_factor<-nrow(de_samples)
height_scale_factor<-as.numeric(nrow(scaled_mat))




h<-Heatmap(scaled_mat,
        
        #cluster_columns = dendsort(hclust(dist(t(scaled_mat)))),
        column_labels = column_labels_mouse_num,
        
        heatmap_legend_param = list(title = ""),
           
        column_title = plot_title,
        
        show_row_names = FALSE,
        show_row_dend = FALSE,
        top_annotation = anno,
        
        
        heatmap_width = unit(.35*width_scale_factor, "in"),
    
        height = unit(.03*height_scale_factor, "in"),
        column_dend_height = unit(.3, "in")

        
        )


gh<-grid.grabExpr(draw(h))


ggsave(paste0("plots/",fs::path_sanitize(paste0("hm-", file_base, "-rpkm.pdf"))),
       plot=gh,
       
       scale = 1.2,
       dpi=600,
       width = unit(.35*width_scale_factor+10, "in"),
       height = unit(.03*height_scale_factor+3, "in"),
       unit="in"
       
       )


print(paste(all_dexps[i], "done"))

}

# stop("end of foreach test loop")

plots<-list.files("plots", full.names = T, pattern="^hm.*\\.pdf$", recursive = F)



qpdf::pdf_combine(plots, output = "plots/combined/hm-combined-1.pdf")





