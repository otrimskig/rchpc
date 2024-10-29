source("libs.R")
library(tidyverse)
library(dtplyr)
suppressPackageStartupMessages(library(ComplexHeatmap))
library(gridtext)


#set list of gene names of interest (mouse format) in dataset.
gene_names<-c("Mdm2", "Cdk4", "Cdkn2a", "Pten", "Atrx", "Nf1")




#read in full rpkm dataset. Filter by specified gene names. 
rpkms<-readRDS("nf1g/ds/vm-02-filtered_rpkms.rds")%>%
  filter(gene_name_ms %in% gene_names)



#pull unique genes found in filtered dataset.
found_genes<-rpkms%>%select(gene_name_ms)%>%unique()%>%pull()


#check that all search terms found a match.
match_check<-tibble("gene_name" = gene_names)%>%
  mutate(matched = if_else(gene_names %in% found_genes, "matched", NA))


#stop progress of code if any unmatched genes detected.
if(nrow(match_check%>%filter(is.na(matched)))==0){
"no unmatched genes. hooray"
}else{
  stop("unmatched genes found")
}
  








#standard imports for nf1 project.



anno_color<-readRDS("nf1g/ds/colors_list.rds")
de_samples<-readRDS("nf1g/ds/v10-per_sample_updated.rds")





#hm annotation setup.

anno<-HeatmapAnnotation(df=de_samples%>%select(#mouse_num, 
                                               patho_grade, 
                                               patho_cat_name,
                                               patho_cat2_name
                                               ), 
                        col=c(anno_color),
                        annotation_name_side = "left",
                        gp = gpar(col = "black"),
                        
                        simple_anno_size = unit(.25, "in"))


#######################################################################

#turn rpkms df to matrix for heatmap. 
mat<-rpkms%>%
  select(gene_name_ms, sample_id, rpkm)%>%
  pivot_wider(names_from = sample_id, values_from = rpkm)%>%
  column_to_rownames("gene_name_ms")%>%
  as.matrix.data.frame()

zmat<-t(scale(t(mat)))

hmz <- Heatmap(zmat,
                
                 top_annotation = anno,
                # column_title = gt_render(
                #   paste0("<span style='font-size:25pt'>Expression Levels Normalized per gene <br>(Relative expression)</span>"), 
                #   r = unit(2, "pt")),
                # 
                # 
                
                
             
                column_title_gp = gpar(font = 2, fontsize = 60),
                row_title_gp = gpar(font = 2, fontsize = 40),
                show_heatmap_legend = FALSE,
                
                
              
                
)

hmz



