source("libs.R")
library(tidyverse)
library(ComplexHeatmap)

###########################################
all_sample_info<-readRDS("k19mf/ds/vm-00-sample_info.rds")

rpkms<-readRDS("k19mf/ds/immune2_rpkms.rds")

stats<-readRDS("k19mf/ds/immune2-acral-v-subq-stats.rds")

groups<-readRDS("k19mf/ds/immune2-genes_groups.rds")
##############################################



#transform from df into a list which is easier to use for subsetting.
groups_list <- groups %>%
  group_by(group) %>%
  summarise(gene_name_ms = list(gene_name_ms)) %>%
  deframe()


group_names<-names(groups_list)


for (g in 1:length(group_names)){


group_name<-group_names[g]

genes_in_group<-groups_list[[group_name]]





int_mat<-rpkms%>%
  
  
  filter(grepl(paste(genes_in_group, collapse = "|"), gene_name_ms))%>%
  
  left_join(all_sample_info%>%select(sample_id,mouse_num))%>%
  select(gene_name_ms, mouse_num, rpkm)%>%
  
    # 
    # 
    # filter(mouse_num!="3126")%>%
  
  
  pivot_wider(names_from = mouse_num, values_from = rpkm)%>%
  column_to_rownames("gene_name_ms")%>%
  as.matrix.data.frame()





int_mat_z<-t(scale(t(int_mat)))




anno<-HeatmapAnnotation(df=all_sample_info%>%  
                          
                          # filter(mouse_num!="3126")%>%
                          
                          select(tumor_type),

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
           
        column_title = group_name,
        column_title_gp = gpar(fontsize = 17, fontface = "bold"),
        
        show_row_names = TRUE,
        show_row_dend = TRUE,
        
        top_annotation = anno,
        
        
        width = unit(.7*ncol(int_mat_z), "cm"),
        height = unit(.5*nrow(int_mat_z), "cm"),
        #column_dend_height = unit(.3, "in")

        
        )


# h1





gh<-grid.grabExpr(draw(h1))


ggsave(paste0("k19mf/plots/", "immune-hm-", group_name, ".pdf"),
       plot=gh,
       
       scale = 1,
       dpi=600,
       
       height=.5*nrow(int_mat_z)+8,
       width =.7*ncol(int_mat_z)+10,
       units = "cm",
       
       limitsize = FALSE
       
       )



}










plots<-fs::dir_info("k19mf/plots", full.names = T, pattern="^hm.*\\.pdf$", recurse = F)%>%tibble()%>%
  filter(modification_time>Sys.time()-lubridate::minutes(5))%>%
  arrange(path)%>%
  filter(type=="file")%>%
  pull(path)



qpdf::pdf_combine(plots, output = "k19mf/plots/combined/immune-hms2.pdf")
