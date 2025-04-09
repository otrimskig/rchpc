source("libs.R")
library(tidyverse)
library(dtplyr)
suppressPackageStartupMessages(library(ComplexHeatmap))
ht_opt$message = FALSE
library(gridtext)
library(circlize)
library(purrr)



vm0<-readRDS("acral_sub_rppa/ds/rppa_df_list0.rds")

sample_info0<-readRDS("acral_sub_rppa/ds/sample_info0.rds")


sample_interest<-sample_info0%>%
  filter(sample_type!="none")

# vm1<-vm0[["L4 (log_2)"]][["df_long"]]

vm1<-vm0[["L4 (linear)"]][["df_long"]]

mat<-vm1%>%
  filter(order %in% sample_interest$order)%>%
  # 
  # filter(order!=455&order!=456)%>%
  # filter(order!=469)%>%
  left_join(sample_interest%>%select(mouse_num, order))%>%
  mutate(sample_name_hm=paste0(mouse_num, "-", order))%>%
  select(-order, -mouse_num)%>%
  pivot_wider(names_from = antibody_name, values_from=rppa_value)%>%
  column_to_rownames("sample_name_hm")%>%
  as.matrix.data.frame()%>%
  t()

mat_num <- apply(mat, 2, as.numeric)


mat_scaled<-t(scale(t(mat_num)))




anno_df<-sample_interest%>%
  select(order, mouse_num, sample_type, a_b)%>%
  mutate(order=paste0("x", order))%>%
  select(-order)
  







anno<-HeatmapAnnotation(df=anno_df, 
                        #col=c(anno_subset),
                        annotation_name_side = "left",
                        
                        gp = gpar(col = "black", fontsize = 3),
                        
                        #annotation_label = c(col_proper_names),
                        annotation_name_gp= gpar(fontsize=7),
                        
                        simple_anno_size = unit(.125, "in"))


# 
# 
#  
# 
#   ha3 = rowAnnotation(a=anno_points(summary_stats_subset[,3],
#                                     pch = c("|"), 
#                                     #gp = gpar(col = 2:3))
#                                     #joyplot_scale = 5,
#                                     width = unit(2, "cm"),
#                                     axis_param = list(
#                                       side = "top")),
#                       show_annotation_name = c(a = FALSE))
#   
  
  

hm1<-Heatmap(mat_scaled,
        top_annotation = anno)

#   
# Heatmap(mat,
#                 #right_annotation = ha3,
#                 #top_annotation = anno,
#                 #column_title = gt_render(
#                 #  paste0("<span style='font-size:25pt'>","Pathway Group: ", subset_factor, "</span><br><span style='font-size:15pt'>expression levels normalized per row (pathway)</span>"), 
#                 #  r = unit(2, "pt")),
#                 
#                 column_dend_height = unit(1, "in"),
#                 row_dend_width = unit(2, "in"),
#                 row_names_side = "right",
#                 row_gap = unit(1, "in"),
#                 
#                 column_title_gp = gpar(font = 2, fontsize = 60),
#                 row_title_gp = gpar(font = 2, fontsize = 40),
#                 show_heatmap_legend = FALSE,
#                 
#                 #width = unit(10, "in"),
#                 #heatmap_height = unit(nrow(gsva_subset)*.3, "in")
#                 
#   )
#   
#   
  

  
#   
#   # 
#   # combo<-draw((hm1), annotation_legend_side = "left")
#   
  gh2<-grid.grabExpr(draw(draw((hm1), annotation_legend_side = "left")))
#   
  ggsave(paste0("acral_sub_rppa/plots/hm-all_samples-l4lin.pdf"),
         plot=gh2,

         scale = 2,
         dpi=600,
         width = 10,
         height = 7,
         unit="in",
         limitsize = FALSE

  )
#   
#   
# }
# 
