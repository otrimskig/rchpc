source("libs.R")
library(tidyverse)
library(dtplyr)
library(ggplot2)
library(ggrepel)
library(ggridges)
library(ggbeeswarm)


dataset<-readRDS("nf1g/atrx_comps/ds/atrx_comps-gsvaNerve sheath tumors.rds")%>%
  rename(fc_orig=fc)%>%
  mutate(fc_abs=abs(fc_orig))%>%
  mutate(log_fc_abs=log(fc_abs))

#comparison set-up
group_a_col<-"coh1_mice_mean_value"
group_b_col<-"coh2_mice_mean_value"



#plot and diff calculator set-up
#IMPORTANT:: assumes z-scores which contain negative values.

save_location<-"nf1g/atrx_comps/plots/" #folder (include trailing slash)


title_of_plot<-"4KO vs. ATRX wt tumors"
subtitle_of_plot<-"Pathway signature comparisons in RNAseq data"


#plot setup.
x_axis_var<-"diff"
y_axis_var<-"pv"
point_labels<-"pn"

x_val_cutoff<-.05 #include log transform.
y_val_cutoff<-(-log10(1E-2)) #include same log transform if nec.




#calculate differences.
dataset2<-dataset%>%
  mutate(diff=!!sym(group_b_col)-!!sym(group_a_col))



dataset3 <- dataset2 %>%
  mutate(
    # Safely apply log transformation with a small value for non-positive values
      # Replace 0 or negative with a small value
    log_pval = -log10(pmax(!!sym(y_axis_var), 1e-10)),  # Ensure p-values are positive for log10
    
    
    significance = case_when(
      abs(!!sym(x_axis_var)) >= abs(x_val_cutoff) & log_pval >= y_val_cutoff ~ "sig_fc",
      
      log_pval >= y_val_cutoff ~ "sig",
      
      abs(!!sym(x_axis_var)) >= abs(x_val_cutoff) ~ "fc",
      
      TRUE ~ "ns"
    ))%>%
  
  mutate(
    label_score_low= if_else(!!sym(x_axis_var)<0, abs(!!sym(x_axis_var))* log_pval, NA),
    label_score_high= if_else(!!sym(x_axis_var)>0, abs(!!sym(x_axis_var))* log_pval, NA),
    
    label_score_rank_low = min_rank(desc(label_score_low)),
    label_score_rank_high = min_rank(desc(label_score_high))
  )%>%
  
  
  mutate(glabel=if_else(significance=="sig_fc"&(label_score_rank_low<=20|label_score_rank_high<=20), !!sym(point_labels), NA))%>%
  
  
  
  mutate(glabel=gsub("_", " ", glabel))%>%

  mutate(glabel=gsub('(.{1,30})(\\s|$)', '\\1\n', glabel))%>%
  mutate(glabel = sub("\n$", "", glabel))%>%
  
  mutate(glabel.a=if_else(!!sym(x_axis_var)<0, glabel, NA))%>%
  mutate(glabel.b=if_else(!!sym(x_axis_var)>0, glabel, NA))
  
  
  
  
  
  

#check max x-value to ensure plot is symmetrical. 
max_abs <- max(abs(dataset3[[sym(x_axis_var)]]), na.rm = TRUE)



# Plot with the updated significance categories
# p<-ggplot(dataset3, aes(x = log_fc, y = log_pval, 
#                      color = significance, label = glabel)) +
  
p<-ggplot(dataset3, aes(x = !!sym(x_axis_var), y = log_pval, 
                          color = significance, label = glabel)) +
  
  
  
  
  
  geom_point(alpha = 0.8, shape=1, size=3) +
  geom_point(alpha = 0.4, shape=19, size=3)+
  # Add vertical lines at fc_cutoff (both directions)
  geom_vline(xintercept = -x_val_cutoff, linetype = "dotted", alpha = 0.7) +
  geom_vline(xintercept = x_val_cutoff, linetype = "dotted", alpha = 0.7) +
  
  # Add a horizontal line at p-value cutoff
  geom_hline(yintercept = y_val_cutoff, linetype = "dotted", alpha = 0.7) +
  
  # Custom color scale for significance categories
  scale_color_manual(values = c("sig_fc" = "purple", 
                                "sig" = "green", 
                                "fc" ="#8aedf2",
                                "ns" = "gray")) +
  
  #expands y-scale to ensure there is room at top of plot. 
  #can remove this if plot has large p-values. 
  scale_y_continuous(expand = expansion(mult = c(0.1, 1)))+
  
  #sets x-axis limits to ensure plot axes limits are symmetrical side-to-side.  
  scale_x_continuous(limits = c(-max_abs, max_abs)) +
  
  geom_label_repel(aes(label=glabel.a),
                  color="black", max.overlaps = 40, 
                  min.segment.length = .1,
                  force=100,
                  xlim  = c(NA, -x_val_cutoff))+
  
  geom_label_repel(aes(label=glabel.b),
                   color="black", max.overlaps = 40, 
                   min.segment.length = .1,
                   force=60,
                   xlim  = c(x_val_cutoff, NA))+
  
  labs(x = "Enrichment Score Difference",
       y = expression("-"~Log[10]~"p-value"),
       title=title_of_plot,
       subtitle = subtitle_of_plot)+
  

  theme_minimal()+
  theme(legend.position = "none",
        plot.title = element_text(size = 25, face = "bold"), 
        plot.subtitle = element_text(size=18),# Increase title size
        axis.title = element_text(size = 16))  # Increase axis label size





ggsave("nf1g/atrx_comps/plots/vol-pathways-atrx_comps_nerve_sheath_tumors.pdf",
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       plot=p,
       limitsize = FALSE,
       
       
       height=15,
       width=25,
       scale = 1,
       dpi=600,
       
       
       
)




