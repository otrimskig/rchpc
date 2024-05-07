setwd("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig")

library(tidyverse)
library(ComplexHeatmap)


heatmap_colors_list<-readRDS("hm_colors_list.rds")


##################################################
#basic inputs.
all_sample_info<-readRDS("23908R/v07-per_sample_info.rds")

ko_rpkms<-readRDS("23908R/v06-all_counts_plus_info.rds")%>%
  filter(gene_name=="Pten"|gene_name== "Cdkn2a"| gene_name=="Nf1"| gene_name=="Atrx"|gene_name=="Nes")

  
ko_genes<-c("Pten", "Cdkn2a", "Nf1","Atrx", "Nes")  

##############################################
#create pre-matrix of counts with all samples and only the 4ko genes. 


mat<-ko_rpkms%>%  
  
  #filter(gene_name==ko_genes[i])%>%
  select(gene_name, mouse_num, rpkm)%>%
  pivot_wider(values_from = rpkm, names_from = mouse_num)%>%
  column_to_rownames("gene_name")

#scale it such that expression is normalize per-row (gene) instead of as a whole. 
scaled_mat<-t(scale(t(mat)))


#create universal annotation for all heatmaps select necessary variables here. 
anno<-HeatmapAnnotation(df=all_sample_info%>%
                          select(aod, resultant_geno, patho_cat, patho_grade),
                        col=c(heatmap_colors_list,
                              aod=circlize::colorRamp2(c(50, 150), c("#001cbd", "white"))),
                        gp = gpar(col = "black")
                        
                        )


# anno<-readRDS("anno.rds")



#saveRDS(anno, "anno.rds")

#make a heatmap with all 4ko genes, with unsupervised clustering.
hm_usc_4ko<-Heatmap(as.matrix(scaled_mat),
            
            column_title = "KO genes with unsupervised clustering",
            show_row_names = TRUE,
            show_row_dend = FALSE,

            top_annotation = anno,
            width = unit(6, "in"), 
            height = unit(2.5, "in")
)

#make a heatmap with all 4ko genes, with unsupervised clustering but separated by genotype.
hm_by_geno_usc_4ko<-Heatmap(as.matrix(scaled_mat),
                    
                    column_title = "KO genes with unsupervised clustering",
                    show_row_names = TRUE,
                    show_row_dend = FALSE,
                    
                    top_annotation = anno,
                    width = unit(6, "in"), 
                    height = unit(2.5, "in")
)




Heatmap(as.matrix(scaled_mat),
        
        column_title = "KO genes with unsupervised clustering",
        show_row_names = TRUE,
        show_row_dend = FALSE,
        
        top_annotation = anno,
        width = unit(6, "in"), 
        height = unit(2.5, "in")
)










stop()





















#make a heat map that is arranged by resultant geno (no clustering)
hm_arr_ge_4ko<-Heatmap(as.matrix(scaled_mat),
            
            column_title = "KO genes arranged by genotype",
            show_row_names = TRUE,

            column_order=all_sample_info%>%arrange(resultant_geno)%>%pull(mouse_num),

            top_annotation = anno,
            
            show_row_dend = FALSE,
            width = unit(6, "in"), 
            height = unit(2.5, "in")
)

hm_arr_aod<-Heatmap(as.matrix(scaled_mat),
                
                column_title = "KO genes arranged by age of death",
                show_row_names = TRUE,
                column_order=all_sample_info%>%arrange(aod)%>%pull(mouse_num),
                
                top_annotation = anno,
                show_row_dend = FALSE,

                width = unit(6, "in"), 
                height = unit(2, "in")
)






#############################################################
#now individual gene heatmaps.






for(i in 1:length(ko_genes)){

mat2<-ko_rpkms%>%  
  
  filter(gene_name==ko_genes[i])%>%
  select(gene_name, mouse_num, rpkm)%>%
  pivot_wider(values_from = rpkm, names_from = mouse_num)%>%
  column_to_rownames("gene_name")


# anno<-HeatmapAnnotation(df=all_sample_info%>%
#                           select(aod, resultant_geno, starts_with("patho")))


husc<-Heatmap(as.matrix(mat2),
              
              column_title = "KO genes: single gene unsupervised clustering",
              show_row_names = TRUE,
              top_annotation = anno,
              
              show_row_dend = FALSE,
              width = unit(6, "in"), 
              height = unit(.5, "in")
)


assign(paste0("hm_unsc_", tolower(ko_genes[i])), husc)
rm(husc)

harr<-Heatmap(as.matrix(mat2),
            
            column_title = "KO genes: single gene arranged by genotype",
            show_row_names = TRUE,

            column_order=all_sample_info%>%arrange(resultant_geno)%>%pull(mouse_num),

            top_annotation = anno,
            show_row_dend = FALSE,
            width = unit(6, "in"), 
            height = unit(.5, "in")
)


assign(paste0("hm_arr_", tolower(ko_genes[i])), harr)
rm(harr)

}


#################################################
#finally, save all to pdfs. 

hm_objects <- ls(pattern = "^hm")
hm_objects

for(o in 1:length(hm_objects)){


gg_out<-grid.grabExpr(draw(get(hm_objects[o])))


ggsave(paste0("plots/", "4ko_", hm_objects[o], "_heatmap.pdf"),
       
       
       plot=gg_out,
       
       scale = 1.2,
       dpi=600,
       width = 10,
       height = 6,
       unit="in"
       
)
}




