source("libs.R")
library(tidyverse)
library(dtplyr)

library(ggplot2)
library(ggrepel)


dataset<-readRDS("nf1g/atrx_comps/ds/atrx_comps-gsvaall_tumor_types.rds")


calculate_signed_fc_tidy <- function(df, group_a_col, group_b_col) {
  df %>%
    mutate(
      signed_fc = ifelse(!!sym(group_a_col) == 0 | !!sym(group_b_col) == 0, 
                         NA, 
                         !!sym(group_b_col) / !!sym(group_a_col)),
      signed_log2_fc = ifelse(!!sym(group_a_col) == 0 | !!sym(group_b_col) == 0, 
                              NA, 
                              sign(!!sym(group_b_col) - !!sym(group_a_col)) * 
                                log2(abs(!!sym(group_b_col)) / abs(!!sym(group_a_col))))
    )
}




dataset<-readRDS("nf1g/atrx_comps/ds/atrx_comps-gsvaall_tumor_types.rds")


dataset<-calculate_signed_fc_tidy(dataset, "coh1_mice_mean_value", "coh2_mice_mean_value")



#dataset<-dataset%>%slice_sample(n = 1000)




title_of_plot<-"4KO vs. ATRX wt tumors"


subtitle_of_plot<-"Pathway signature comparisons in RNAseq data"



save_location<-"nf1g/atrx_comps/plots/" #folder (include trailing slash)

x_axis_var<-"fc"

y_axis_var<-"pv"

point_labels<-"pn"


fc_cutoff<-1.25

fc_cutoff<-1.5

p_value_cutoff<-1E-2





dataset2 <- dataset %>%
  mutate(
    # Safely apply log transformation with a small value for non-positive values
    log_fc = log(pmax(!!sym(x_axis_var), 1e-10)),  # Replace 0 or negative with a small value
    log_pval = -log10(pmax(!!sym(y_axis_var), 1e-10)),  # Ensure p-values are positive for log10
    significance = case_when(

      abs(log_fc) >= abs(log(fc_cutoff)) & log_pval >= -log10(p_value_cutoff) ~ "sig_hfc",
      log_pval >= -log10(p_value_cutoff) ~ "sig",
      TRUE ~ "ns"
    ))%>%
      mutate(glabel=if_else(significance=="sig_hfc"|significance=="sig", !!sym(point_labels), NA))%>%

      abs(log_fc) >= abs(log(fc_cutoff)) & log_pval >= -log10(p_value_cutoff) ~ "sig_fc",
      log_pval >= -log10(p_value_cutoff) ~ "sig",
      abs(log_fc) >= abs(log(fc_cutoff)) ~ "fc",
      TRUE ~ "ns"
    ))%>%
      mutate(glabel=if_else(significance=="sig_fc"|significance=="sig", !!sym(point_labels), NA))%>%

      mutate(glabel=gsub("_", " ", glabel))%>%
  #mutate(glabel = str_replace_all(glabel, "(.{1,25})(?=\\s)", "\\1\n"))
  mutate(glabel=gsub('(.{1,30})(\\s|$)', '\\1\n', glabel))%>%
  mutate(glabel = sub("\n$", "", glabel))
  


#check max x-value to ensure plot is symmetrical. 
max_abs <- max(abs(log(dataset2[[x_axis_var]])), na.rm = TRUE)



# Plot with the updated significance categories
p<-ggplot(dataset2, aes(x = log(!!sym(x_axis_var)), y = -log10(!!sym(y_axis_var)), 
                     color = significance, label = glabel)) +
  geom_point(alpha = 0.8, shape=1, size=3) +
  geom_point(alpha = 0.4, shape=19, size=3) +
  # Add vertical lines at fc_cutoff (both directions)
  geom_vline(xintercept = log(fc_cutoff), linetype = "dotted", alpha = 0.7) +
  geom_vline(xintercept = -log(fc_cutoff), linetype = "dotted", alpha = 0.7) +
  
  # Add a horizontal line at p-value cutoff
  geom_hline(yintercept = -log10(p_value_cutoff), linetype = "dotted", alpha = 0.7) +
  
  # Custom color scale for significance categories

  scale_color_manual(values = c("sig_hfc" = "purple", 
                                "sig" = "green", 
                                "ns" = "gray")) +
  
  
  
  geom_label_repel(color="black", max.overlaps = 50, 
                  min.segment.length = .1,
                  force=20)+
  
  labs(x = expression(Log[2]~"Fold Change"),
       y = expression("-"~Log[10]~"p-value"))+
  
  ggtitle(title_of_plot)+
  theme_minimal()+
  theme(legend.position = "none")

  scale_color_manual(values = c("sig_fc" = "purple", 
                                "sig" = "green", 
                                "fc" ="#8aedf2",
                                "ns" = "gray")) +
  
  #expands y-scale to ensure there is room at top of plot. 
  #can remove this if plot has large p-values. 
  scale_y_continuous(expand = expansion(mult = c(0.1, 1)))+
  
  #sets x-axis limits to ensure plot axes limits are symmetrical side-to-side.  
  scale_x_continuous(limits = c(-max_abs, max_abs)) +
  
  geom_label_repel(color="black", max.overlaps = 50, 
                  min.segment.length = .1,
                  force=10)+
  
  labs(x = expression(Log[2]~"Fold Change"),
       y = expression("-"~Log[10]~"p-value"),
       title=title_of_plot,
       subtitle = subtitle_of_plot)+
  

  theme_minimal()+
  theme(legend.position = "none",
        plot.title = element_text(size = 25, face = "bold"), 
        plot.subtitle = element_text(size=18),# Increase title size
        axis.title = element_text(size = 16))  # Increase axis label size






ggsave("nf1g/atrx_comps/plots/vol-pathways-atrx_comps_all_tumor_types.pdf",
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       plot=p,
       limitsize = FALSE,
       
       
       height=12,
       width=25,
       scale = 1,
       dpi=600,
       
       
       
)





stop()

EnhancedVolcano(dataset,
                lab = dataset$pn,
                x = 'fc',
                y = 'pv',
                title = 'N061011 versus N61311',
                pCutoff = 10e-32,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)

