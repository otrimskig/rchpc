source("libs.R")
library(tidyverse)
library(dtplyr)
suppressPackageStartupMessages(library(ComplexHeatmap))
library(gridtext)


#set list of gene names of interest (mouse format) in dataset.
gene_names<-c("Trp53", "Rb1", "Mdm2", "Cdk4", "Cdkn2a", "Cdkn2b", "Pten", "Atrx", "Nf1", "Fat1", "Hdlbp")




#read in full rpkm dataset. Filter by specified gene names. 
rpkms<-readRDS("nf1g/ds/vm-02-filtered_rpkms.rds")%>%
  mutate(gene_name_lower=tolower(gene_name_ms))%>%
  filter(gene_name_lower %in% tolower(gene_names))



#pull unique genes found in filtered dataset.
found_genes<-rpkms%>%select(gene_name_ms)%>%unique()%>%pull()


#check that all search terms found a match.
match_check<-tibble("gene_name" = gene_names)%>%
  mutate(matched = if_else(gene_names %in% found_genes, "matched", NA))


#stop progress of code if any unmatched genes detected.
if(nrow(match_check%>%filter(is.na(matched)))==0){
"no unmatched genes. hooray"
}else{
  
  no_matches<-match_check%>%
    filter(is.na(matched))%>%
    pull(gene_name)
  
  stop(paste0("unmatched genes found: ", no_matches))
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
                
               height=unit(1*nrow(zmat),"cm"),
             
                column_title_gp = gpar(font = 2, fontsize = 60),
                row_title_gp = gpar(font = 2, fontsize = 40),
                show_heatmap_legend = FALSE
                
                
              
                
)

# hmz



gh<-grid.grabExpr(draw(hmz))


# ggsave(paste0("nf1g/plots/", "hm-single-genes.pdf"),
#        plot=gh,
#        
#        scale = 1,
#        dpi=600,
#        
#        height=.5*nrow(zmat)+8,
#        width =.7*ncol(zmat)+10,
#        units = "cm",
#        
#        limitsize = FALSE
#        
# )




hm_list<-c("Combined" = hmz)


single_genes<-rownames(zmat)




for (sg in 1:length(single_genes)){

  
sgn<-single_genes[sg]  


zmat_single<-zmat[sgn, ,drop = FALSE]


hm_list[[sgn]] <- Heatmap(zmat_single,
               
               top_annotation = anno,
               # column_title = gt_render(
               #   paste0("<span style='font-size:25pt'>Expression Levels Normalized per gene <br>(Relative expression)</span>"), 
               #   r = unit(2, "pt")),
               # 
               # 
               height=unit(1,"cm"),
               
               
               column_title_gp = gpar(font = 2, fontsize = 60),
               row_title_gp = gpar(font = 2, fontsize = 40),
               show_heatmap_legend = FALSE
               
               
               
               
)

  
}








library(ComplexHeatmap)
library(grid)
library(cowplot)

# Assuming hm_list is your list of ComplexHeatmap objects
# Capture each heatmap in hm_list as a grob
plot_grid_list <- lapply(hm_list, function(hm) {
  grid.grabExpr(draw(hm))
})


final_plot <- plot_grid(
  plotlist = plot_grid_list, 
  ncol = 1, 
  align = "v",
  rel_heights = rep(1, length(plot_grid_list))  # Adjust this vector to fine-tune spacing
)


output_loc<-paste0("nf1g/plots/", "hm-single-genes.pdf")

ggsave(output_loc,
       plot=final_plot,
       height=100,
       width=30,
       scale = 1,
       dpi=600,
       limitsize = FALSE
)

Cairo::CairoPDF(file = output_loc, 
         title = "Selected genes from NF1 glioma RNA seq",
         height=100,
         width=30)
print(final_plot)
dev.off()
