setwd("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig")

library(tidyverse)
library(ComplexHeatmap)
library(EnhancedVolcano)

###########################################

logfc_threshold<-2
fdr_threshold<-.001

rds_file_path<-"dexps/dexp-patho_cat-NEDv3.rds"





##############################################


plot_title<-basename(rds_file_path)%>%
  sub(".rds$", "", .)%>%
  sub("^dexp-", "", .)%>%
  sub("-", ": ", .)%>%
  sub("v", " vs. ", .)



file_base<-basename(rds_file_path)%>%
  sub(".rds$", "", .)%>%
  sub("^dexp-", "", .)




de_df<-readRDS(rds_file_path)


all_sample_info<-readRDS("23908R/v07-per_sample_info.rds")

#filter based on only samples included in DE groups.
de_samples<-de_df%>%
  colnames()%>%
  as_tibble()%>%
  rename(sample_id=value)%>%
  filter(grepl("^rpkm",sample_id))%>%
  mutate(sample_id=sub("^rpkm_","", sample_id))%>%
  
  left_join(all_sample_info)%>%
  mutate(patho_cat=as_factor(patho_cat))%>%
  mutate(patho_cat2=as_factor(patho_cat2))



hm_mat<-de_df%>%
  rename_with(~ sub("^rpkm_", "", .), starts_with("rpkm_"))%>%
 
  filter(abs(logFC)>=logfc_threshold)%>%
  filter(FDR<fdr_threshold)%>%
  
  select(-c(gene_id:FDR))%>%
  column_to_rownames("gene_name")



if (nrow(hm_mat) > 500) {
  # If true, sample 500 rows
  hm_mat_500 <- hm_mat %>% sample_n(500)
} else {
  # If false, don't sample
  hm_mat_500 <- hm_mat
}


scaled_mat<-t(scale(t(hm_mat_500)))





anno<-HeatmapAnnotation(df=de_samples%>%select(patho_grade, patho_cat, patho_cat2, patho_cat_det,aod), 
                        # col= list(resultant_geno=c("nf1 KO; pten KO; ink KO; atrx KO"="#66c2a5", 
                        #                            "nf1 KO; pten KO; ink KO; atrx wt"="#fc8d62", 
                        #                            "nf1 wt; pten KO; ink KO; atrx KO"="navy", 
                        #                            "nf1 wt; pten wt; ink wt; atrx wt"="#e7298a")
                        
                        gp = gpar(col = "black"))



width_scale_factor<-nrow(de_samples)
height_scale_factor<-nrow(scaled_mat)

h<-Heatmap(scaled_mat,
        
        column_title = plot_title,
        
        show_row_names = FALSE,
        show_row_dend = FALSE,
        top_annotation = anno,
        
        
        heatmap_width = unit(.35*width_scale_factor, "in"),
        heatmap_height = unit(.02*height_scale_factor, "in")
        
        
        )



gh<-grid.grabExpr(draw(h))


ggsave(paste0("plots/",fs::path_sanitize(paste0(file_base, "-rpkm-heatmap.pdf"))),
       plot=gh,
       
       scale = 1.2,
       dpi=600,
       width = 10,
       height = 12,
       unit="in"
       
       )

# plot(gh)




# gb = grid.grabExpr(draw(h))
# gb
# is.grob(gb)
# 
# 
# 
# pushViewport(viewport(width = 1, height = 1))
# grid.draw(gb)




########################################
#volcano plot for same data.



vdf<-de_df

p<-EnhancedVolcano(vdf, 
                lab=vdf$gene_name,
                x="logFC",
                y="PValue",
                axisLabSize = 12,
                labSize = 5.0,
                drawConnectors = TRUE,
                boxedLabels = FALSE,
                colConnectors = "grey50",
                arrowheads = FALSE,
                min.segment.length=.3,
                title=plot_title,
                subtitle = "",
                maxoverlaps = 10,
                #typeConnectors="open",
                
                endsConnectors="first",
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 4.0
                
                
                #selectLab = c("Pten", "Cdkn2a", "Nf1", "Atrx")
                
                )


ggsave(paste0("plots/",fs::path_sanitize(paste0(file_base, "-volcano.pdf"))),
       plot=p,
       
       scale = 1.2,
       dpi=600,
       width = 10,
       height = 8,
       unit="in"
       
)


