source("libs.R")
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(ggbeeswarm)

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

###########################################
sample_info<-readRDS("nf1g/ds/v10-per_sample_updated.rds")

gsva_pathway_stats<-readRDS("nf1g/gsvas/gsva_top500_output.rds")$results2%>%
  arrange(desc(abs_diff))%>%
  slice(1:100)

gsva_u0<-readRDS("nf1g/ds/gsva/gsva_pathways_matrix.rds")[gsva_pathway_stats$Pathway,]



samples<-sample_info%>%
  filter(patho_cat_name=="Gliomas")%>%
  pull(sample_id)


gsva_u1<-gsva_u0[,samples]


gsva_z1<-t(scale(t(gsva_u1)))






timex_ds<-gsva_u1%>%
  as_tibble(., rownames = "gene_set")%>%
  pivot_longer(starts_with("x"), names_to="sample_id", values_to = "enrichment_score")%>%
  left_join(sample_info)%>%

  group_by(gene_set, resultant_geno)%>%
  mutate(outlier = ifelse(is_outlier(enrichment_score), enrichment_score, as.numeric(NA)))%>%
  ungroup()



p <- timex_ds %>% 
  arrange(gene_set) %>%
  #slice(1:108) %>%
  ggplot(aes(x = resultant_geno, y = enrichment_score, color = patho_grade)) + 
  facet_wrap(vars(gene_set)) +
  geom_beeswarm(size = 1, alpha = .7, cex = 6, stroke = 1.5) +  # Removed color = "black"
  theme_classic() +
  scale_color_manual(values = c("2" = "#7fc97f", "3" = "#fdc086", "4" = "#beaed4"))+
  
  
  
  
  theme(
    #legend.position = "none",  # Hide legend
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    #plot.margin = margin(l = 200, unit = "pt") # Adjust aspect ratio if needed
  )









print(p)
  

  facet_wrap(~gene_set, scales = "free_x", nrow=3) +

  
  scale_x_discrete(expand = expansion(add = c(2, 2)), labels = c("", ""))+
  

  xlab("")+

  scale_fill_manual(values=c("subq"="red", "acral"= "blue"))+

  
  geom_text(aes(label=ifelse(sample_type=="acral"&!is.na(outlier),
                             as.character(mouse_num),'')),hjust=1.5,vjust=0.5,
            size = 3)+
  
  geom_text(aes(label=ifelse(sample_type=="subq"&!is.na(outlier),
                             as.character(mouse_num),'')),hjust=-.75,vjust=0.5,
            size = 3)+
  
  
  
  
  
  
  theme_classic() +
  theme(
    #legend.position = "none",  # Hide legend
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    #plot.margin = margin(l = 200, unit = "pt") # Adjust aspect ratio if needed
  ) +
  #coord_cartesian(clip = "off")+
  
  labs(
    x = "Gene Set",
    y = "Normalized Enrichment Score",
    color = "Sample Type"
  )




















stop()


ggsave("acral_paired/plots/cell_types_matched_scatter_labeled_normalized.pdf",
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       plot=p,
       limitsize = FALSE,
       
       
       height=10,
       width=20,
       scale = 2,
       dpi=600,
       
       
       
)










