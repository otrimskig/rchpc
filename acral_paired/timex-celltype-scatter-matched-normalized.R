source("libs.R")
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)



###########################################
sample_info<-readRDS("acral_paired/ds/v00-sample_info.rds")


gsva_ds<-readRDS("acral_paired/ds/gsva_z.rds")



is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}






timex_ds<-gsva_ds%>%
  as_tibble(., rownames = "gene_set")%>%
  filter(grepl("^timex_", gene_set))%>%
  mutate(gene_set=gsub("^timex_", "", gene_set))%>%
  pivot_longer(ends_with("subq")|ends_with("acral"), names_to="mouse_num_type_clean", values_to = "enrichment_score")%>%
  
  left_join(readRDS("acral_paired/ds/name_map.rds"))%>%
  group_by(sample_id)%>%
  ungroup()%>%
  arrange(sample_type)%>%
  group_by(gene_set, sample_type)%>%


  mutate(outlier = ifelse(is_outlier(enrichment_score), enrichment_score, as.numeric(NA)))%>%
  ungroup()




gene_set_index<-timex_ds%>%
  group_by(gene_set)%>%
  slice(1)%>%
  ungroup()%>%
  mutate(gene_set_index=1:n())







# timex_ds2<-timex_ds%>%
# 
# 
# 
# 
# timex_ds$sample_type_factor <- factor(timex_ds$sample_type, levels = c("subq", "acral"))
p<-timex_ds %>% 
  mutate(new = factor(interaction(gene_set, sample_type))) %>%
  
  
  ggplot(aes(x = new, y = enrichment_score)) + 
  #geom_boxplot(outlier.shape = NA, width=.4) + 
  
  geom_line(aes(group = mouse_num), color = "gray") +
  geom_dotplot(aes(fill = sample_type),
               alpha=.7,
               binaxis = "y", 
               #binwidth = 1E-3, 
               stackdir = "center", 
               stackratio = .5,
               dotsize=1,
               position = position_dodge(20)) +
 
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










