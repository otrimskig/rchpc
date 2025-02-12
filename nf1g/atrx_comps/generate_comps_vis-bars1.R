source("libs.R")
library(tidyverse)
library(dtplyr)
library(ggplot2)
library(ggrepel)
library(ggridges)
library(ggbeeswarm)


dataset<-readRDS("nf1g/atrx_comps/ds/atrx_comps-gsvaSpindle and epithelioid tumors.rds")%>%
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
bar_label<-"pn"

#x_val_cutoff<-.05 #include log transform.
p_val_cutoff<-(-log10(1E-2)) #include same log transform if nec.




#calculate differences.
dataset2<-dataset%>%
  mutate(diff=!!sym(group_b_col)-!!sym(group_a_col))



dataset3 <- dataset2 %>%
  
  mutate(glabel=gsub("_", " ", !!sym(bar_label)))%>%
  
  #mutate(glabel=if_else(nchar(glabel)>30, paste0(substr(glabel, 1,30), "..."), glabel))%>%
  
  
  filter(!grepl("^HP ", glabel))%>%
  filter(!grepl("^MIR[0-9]", glabel))%>%

  
  mutate(pval_rank= min_rank((pv)))%>%
  
  
  filter(pval_rank<=50)%>%
  arrange(desc(diff))%>%
  mutate(glabel = factor(glabel, levels = rev(glabel)))%>%
  mutate(sign=if_else(diff>0, "1","0"))
  


ggplot(dataset3, aes(x = diff, y = glabel, fill=sign)) +
  geom_col(width=.05)+
  geom_vline(xintercept = 0, alpha = 0.7)+
  theme_classic()
  
  
  
  
  
  
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





ggsave("nf1g/atrx_comps/plots/wf-pathways-atrx_comps_all_tumors.pdf",
       
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




