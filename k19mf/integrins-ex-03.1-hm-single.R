source("libs.R")
library(tidyverse)
library(ComplexHeatmap)
# library(EnhancedVolcano)
# library(circlize)


###########################################
#set threshold values for deexps datasets.

logfc_threshold<-2
fdr_threshold<-.001

all_sample_info<-readRDS("k19mf/ds/vm-00-sample_info.rds")
reads<-readRDS("k19mf/ds/intds.rds")

##############################################




int_mat<-reads%>%
  left_join(all_sample_info%>%select(sample_id,mouse_num))%>%
  select(gene_name_ms, mouse_num, rpkm)%>%
  
  

    filter(mouse_num!="3126")%>%
  
  
  pivot_wider(names_from = mouse_num, values_from = rpkm)%>%
  column_to_rownames("gene_name_ms")%>%
  as.matrix.data.frame()










int_mat_z<-t(scale(t(int_mat)))





anno<-HeatmapAnnotation(df=all_sample_info%>%  
                          filter(mouse_num!="3126")%>%
                          select(tumor_type)%>%mutate(tumor_type=gsub("foot pad", "acral", tumor_type)),

                        col=list(tumor_type=c("subq" = "#fffa99", "acral"="#99fff6")),

                        show_annotation_name = FALSE,
                        annotation_legend_param = list(tumor_type = list(title = "Injection Site")),
                        
                        
                        gp = gpar(col = "black"),

                        simple_anno_size = unit(.125, "in"))


# 
# width_scale_factor<-nrow(de_samples)
# height_scale_factor<-as.numeric(nrow(scaled_mat))
# 
# height_scale_factor




h1<-Heatmap(int_mat_z,
        
        #cluster_columns = dendsort(hclust(dist(t(scaled_mat)))),
        #column_labels = column_labels_mouse_num,
        
        heatmap_legend_param = list(title = "z-score"),
           
        #column_title = plot_title,
        
        show_row_names = TRUE,
        show_row_dend = TRUE,
        
        top_annotation = anno
        
        
        #heatmap_width = unit(.35*width_scale_factor, "in"),
    
        #height = unit(.03*height_scale_factor, "in"),
        #column_dend_height = unit(.3, "in")

        
        )


h1








gh<-grid.grabExpr(draw(h1))


ggsave(paste0("k19mf/plots/", "integrins-hm-ex.pdf"),
       plot=gh,
       
       scale = 2,
       dpi=600
       
       )








