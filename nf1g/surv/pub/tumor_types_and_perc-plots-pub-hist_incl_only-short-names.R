source("libs.R")
suppressMessages(source("nf1g/colors_map_create.R", echo = FALSE))

library(tidyverse)
library(ggpubr)
library(gridExtra)
library(ggkm)
library(scales)
library(survival)
library(survminer)
library(RColorBrewer)

source("ggplot_draw_square.R")


########dataset read in and construction######################
#reading in current dataset. 
coh1<-readRDS("nf1g/surv/cohorts-2025-03-07-v2.rds")
col_map<-readRDS("nf1g/surv/colors_map_surv.rds")



df1<-coh1%>%
#most conservative exclusion criteria
filter(is.na(exclude))%>%
  filter(!is.na(include_in_surv))%>%
  filter(is.na(exclude_hist))

aspectratio<-.6


# genotypes<-df1%>%select(resultant_geno)%>%unique()%>%arrange()%>%pull()



df_props<-df1%>%

  select(resultant_geno, hist_cat_name)%>%
  group_by(resultant_geno,hist_cat_name)%>%
  
  #count occurrences of each hist_cat per geno.
  reframe(n=n() )%>%
  
  #calculate percentage/likelihood of hist_cat per geno.
  group_by(resultant_geno)%>%
  reframe(hist_cat_name, n=n,total_n=sum(n))%>%
  mutate(perc=n/total_n*100)%>%

  #add in a dummy number for each hist_cat-geno combo to ensure consistency in plot display groups.
  complete(resultant_geno, hist_cat_name, fill = list(perc = -.5))%>%
  
  #get back info for resultant_geno_name ("proper" cohort names), to use if desired.
  left_join(coh1%>%select(resultant_geno, resultant_geno_name)%>%unique())%>%
  mutate(resultant_geno = factor(resultant_geno, levels = names(col_map$resultant_geno)))%>%
  mutate(resultant_geno_name = factor(resultant_geno_name, levels = names(col_map$resultant_geno_name)))%>%
  mutate(hist_cat_name = factor(hist_cat_name, levels = names(col_map$hist_cat_name)))



df_props$hist_cat_name_numeric <- as.numeric(as.factor(df_props$hist_cat_name))*0.6
#cat_levels <- unique(df_props$hist_cat_name)
#df_props$hist_cat_name_numeric <- match(df_props$hist_cat_name, cat_levels) * 0.6  # Adjust 0.8 to fine-tune spacing



p2<-ggplot(df_props) +
  geom_col(aes(x = hist_cat_name_numeric, 
               fill = resultant_geno,
               y = perc),
           width = 0.4,  
           position = position_dodge(0.5),
           key_glyph = draw_square)  +
  theme_classic()+




  scale_fill_manual(values=col_map$resultant_geno,
                    
                    
                   labels = scales::label_wrap(20))+ #wrap long strings into multiple lines.
  
  theme(
    axis.text.x = element_text(size=12,angle = 45, hjust = 1),
    plot.margin = margin(100, 100, 100, 100))+
  ggtitle("Cohort prevalence by tumor type")+
  labs(fill=NULL,
       x=NULL,
       y="% of Each Cohort")+
  
  
  
  scale_x_continuous(breaks = unique(df_props$hist_cat_name_numeric),
                     labels = unique(df_props$hist_cat_name))+
  theme(
    legend.text = element_text(size = 8, hjust = 0, vjust=0.5), 
    #legend.key.height = unit(5, "mm"),
    legend.key.spacing.y = unit(5, 'mm'),
    aspect.ratio=.6)+
  
guides(fill = guide_legend(byrow = TRUE))
  
  
#p2




metadata_text <- paste0("src: ", 
                        rstudioapi::getSourceEditorContext()$path %>%
                          sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.) %>%
                          sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.), 
                        "\n", "time: ", 
                        lubridate::round_date(Sys.time(), "second"))

text_grob <- textGrob(metadata_text, x=.05, just="left", gp = gpar(fontsize = 9, col = "gray30"))

p2src<-grid.arrange(p2, text_grob, ncol = 1, heights = c(3, 0.3))

#p2src


ggsave("nf1g/surv/pub/tumor_types-perc-noXH-v0-geno.pdf",
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       plot=p2src,
       limitsize = FALSE,
       
       
       height=10,
       width=15,
       scale = 1,
       dpi=600,
       
       
       
)  






df_props$resultant_geno_name_numeric <- as.numeric(as.factor(df_props$resultant_geno_name))
geno_levels <- unique(df_props$resultant_geno_name)
df_props$resultant_geno_name_numeric <- match(df_props$resultant_geno_name, geno_levels) * 0.6  # Adjust 0.8 to fine-tune spacing






p3<-ggplot(df_props) +
  geom_col(aes(x = resultant_geno_name_numeric, 
               fill = hist_cat_name,
               y = perc),
           width = 0.4,  
           position = position_dodge(0.5)) +  
  theme_classic()+
  
  
  
  
  scale_fill_manual(values=col_map$hist_cat_name)+
  
  theme(
    axis.text.x = element_text(size=12,angle = 45, hjust = 1),
    plot.margin = margin(100, 100, 100, 100),
    aspect.ratio=.6
  ) +
  ggtitle("Tumor type prevalence by cohort")+
  labs(fill=NULL,
       x=NULL,
       y="% tumor incidence")+
  
  
  
  scale_x_continuous(
    breaks = unique(df_props$resultant_geno_name_numeric),
    labels = str_wrap(unique(df_props$resultant_geno), width = 18)  # Wrap long labels
  )+
  
  
  geom_segment(aes(
    x = as.numeric(resultant_geno_name_numeric) - 0.1,   # Slightly offset from bars
    xend = as.numeric(resultant_geno_name_numeric) + 0.1, # Slightly offset from bars
    y = -2.5,   # Adjust as needed for the vertical position of the line
    yend = -2.5,  # Keep line horizontal at the same position
    color = resultant_geno  # Color the line based on resultant_geno_name
  ), linewidth = 1.5) +  # Line width
  scale_color_manual(values = col_map$resultant_geno) +
  
  guides(color = "none") 



metadata_text <- paste0("src: ", 
                        rstudioapi::getSourceEditorContext()$path %>%
                          sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.) %>%
                          sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.), 
                        "\n", "time: ", 
                        lubridate::round_date(Sys.time(), "second"))

text_grob <- textGrob(metadata_text, x=.05, just="left", gp = gpar(fontsize = 9, col = "gray30"))

p3src<-grid.arrange(p3, text_grob, ncol = 1, heights = c(3, 0.3))






ggsave("nf1g/surv/pub/tumor_types-perc2-noXH-v0-geno.pdf",
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       plot=p3src,
       limitsize = FALSE,
       
       
       height=10,
       width=15,
       scale = 1,
       dpi=600,
       
       
       
)  

