source("libs.R")

library(ggplot2)
library(ggrepel)
library(ggridges)
library(ggbeeswarm)


rds_files<-list.files("nf1g/custom_diff/dexps", full.names = T)


save_location<-"nf1g/custom_diff/plots/" #folder (include trailing slash)
output_file_suffix<-"label_all_mans02"



x_val_cutoff<-.5 #include log transform.


use_FDR_cutoff<-.05 #set to NULL if using a manual cutoff
y_val_cutoff<-NULL #include same log transform if nec.


x_axis_var<-"logFC"
#y_axis_var<-"PValue"
point_labels<-"gene_name_ms"



manual_additions<-list(starts_with=c("Parp", "rad", "idh", "ifit",
                                    "p2r", "selp", "cx3cr", "cd", "tag"), 
                       contains=c("mre11"),
                       ends_with=c(NA),
                       full_gene_names=c(NA))

keep_only_sigfc_manual_highlights<-T





for(rd in 1:length(rds_files)){
cat("\n\nworking on ", cyan(rds_files[rd]), "\n")

diff0<-readRDS(rds_files[rd])
#extract main dataset from list
diff1<-diff0$data

title_of_plot<-paste0(diff0$meta$comparison$groupa, " vs. ", diff0$meta$comparison$groupb)
subtitle_of_plot<-paste0(diff0$meta$subset1_name)



if(is.null(use_FDR_cutoff)&is.null(y_val_cutoff)){
  
  cat(red("use_FDR_cutoff and y_val_cutoff cannot both be null.\n"))
  stop()
}

if(!is.null(use_FDR_cutoff)){
  
  #this selects the value closest to your set FDR.
  y_val_cutoff<-diff1 %>%
    slice_min(order_by = abs(FDR - use_FDR_cutoff), n = 1)%>%
    slice(1)%>%
    pull(PValue)%>% #get corresponding log10 pval to use for cutoff
    -log10(.)%>%   
     round(digits = 2) #round to a reasonable amt of digits.
  
}



#transform to get labels and colors by sig.
diff2<-diff1%>%
  mutate(
    # Safely apply log transformation with a small value for non-positive values
    # Replace 0 or negative with a small value
    neg_log10_pval = -log10(pmax(PValue, 1e-10)))%>%# Ensure p-values are positive for log10
  relocate(neg_log10_pval, .after="logFC")%>%
  
    
  mutate(significance = case_when(
      abs(!!sym(x_axis_var)) >= abs(x_val_cutoff) & neg_log10_pval >= y_val_cutoff ~ "sig_fc",
      
      neg_log10_pval >= y_val_cutoff ~ "sig",
      
      abs(!!sym(x_axis_var)) >= abs(x_val_cutoff) ~ "fc",
      
      TRUE ~ "ns"
    ))%>%
  
  
 
  mutate(
    manual_highlights = if_else(
      (!is.na(manual_additions$starts_with) & str_starts(gene_name_ms, regex(str_c(manual_additions$starts_with, collapse = "|"), ignore_case = TRUE))) |
        (!is.na(manual_additions$contains) & str_detect(gene_name_ms, regex(str_c(manual_additions$contains, collapse = "|"), ignore_case = TRUE))) |
        (!is.na(manual_additions$ends_with) & str_ends(gene_name_ms, regex(str_c(manual_additions$ends_with, collapse = "|"), ignore_case = TRUE))) |
        (!is.na(manual_additions$full_gene_names) & gene_name_ms %in% manual_additions$full_gene_names),
      gene_name_ms,
      NA_character_
    )
  )





if (keep_only_sigfc_manual_highlights) {
  diff3 <- diff2 %>%
    mutate(manual_highlights = if_else(significance == "sig_fc", manual_highlights, NA))
} else {
  diff3 <- diff2
}
  
  
  
diff4<-diff3%>%


  mutate(significance=if_else(!is.na(manual_highlights), paste0(significance, "--manh"), significance))%>%
  
  
  mutate(
    label_score_low= if_else(!!sym(x_axis_var)<0, abs(!!sym(x_axis_var))* neg_log10_pval, NA),
    label_score_high= if_else(!!sym(x_axis_var)>0, abs(!!sym(x_axis_var))* neg_log10_pval, NA),
    
    label_score_rank_low = min_rank(desc(label_score_low)),
    label_score_rank_high = min_rank(desc(label_score_high))
  )%>%
  
  
  mutate(glabel=if_else(significance=="sig_fc"&(label_score_rank_low<=20|label_score_rank_high<=20), !!sym(point_labels), NA))%>%
  mutate(glabel=if_else(!is.na(manual_highlights), gene_name_ms, glabel))%>%
 
  
  mutate(glabel=gsub("_", " ", glabel))%>%
  
  # mutate(glabel=gsub('(.{1,30})(\\s|$)', '\\1\n', glabel))%>%
  # mutate(glabel = sub("\n$", "", glabel))%>%
  # 
  mutate(glabel.a=if_else(!!sym(x_axis_var)<0, glabel, NA))%>%
  mutate(glabel.b=if_else(!!sym(x_axis_var)>0, glabel, NA))
  
  
#check max x-value to ensure plot is symmetrical. 
max_abs <- max(abs(diff4[[sym(x_axis_var)]]), na.rm = TRUE)




cat("generating plot...\n")

vol1<-ggplot(diff4, aes(x = logFC, y = neg_log10_pval, 
                     color = significance, label = glabel)) +
  
  
  
  
  
  geom_point(data = diff4%>%filter(is.na(manual_highlights))%>%
                                     filter(significance!="sig_fc"),
             alpha = 0.8, shape=1, size=3) +
  geom_point(data = diff4%>%filter(is.na(manual_highlights))%>%
                                     filter(significance!="sig_fc"),
             alpha = 0.4, shape=19, size=3)+
  
  
  geom_point(data = diff4%>%filter(is.na(manual_highlights))%>%
                                     filter(significance=="sig_fc"),
             alpha = 0.8, shape=1, size=3) +
  geom_point(data = diff4%>%filter(is.na(manual_highlights))%>%
                                     filter(significance=="sig_fc"),
             alpha = 0.4, shape=19, size=3)+
  
  
 
  
  geom_point(data = diff4%>%filter(!is.na(manual_highlights)),
             alpha = .7, shape=19, size=3,
             color="#ff7f50") +
  geom_point(data = diff4%>%filter(!is.na(manual_highlights)),
             alpha = 1, shape=1, size=3,
             color="black") +

  
  
  
  


  # Add vertical lines at fc_cutoff (both directions)
  geom_vline(xintercept = -x_val_cutoff, linetype = "dotted", alpha = 0.7) +
  geom_vline(xintercept = x_val_cutoff, linetype = "dotted", alpha = 0.7) +
  
  # Add a horizontal line at p-value cutoff
  geom_hline(yintercept = y_val_cutoff, linetype = "dotted", alpha = 0.7)+
  
  # Custom color scale for significance categories
  scale_color_manual(values = c("ns" = "gray",
                                "fc" ="#8aedf2",
                                "sig" = "green", 
                          
                                "sig_fc" = "purple"
                               
                                
                                ))+
  
 
  
  
  
  #expands y-scale to ensure there is room at top of plot. 
  #can remove this if plot has large p-values. 
  scale_y_continuous(expand = expansion(mult = c(0.1, 1)))+
  
  #sets x-axis limits to ensure plot axes limits are symmetrical side-to-side.  
  scale_x_continuous(limits = c(-max_abs, max_abs)) +
  
  geom_label_repel(data = diff4%>%filter(is.na(manual_highlights)),
                   aes(label=glabel.a),
                   color="black", max.overlaps = 15, 
                   min.segment.length = .1,
                   point.padding=.5,
                   force=5,
                   xlim  = c(NA, -x_val_cutoff))+
  
  geom_label_repel(data = diff4%>%filter(is.na(manual_highlights)),
                   aes(label=glabel.b),
                   color="black", max.overlaps = 15, 
                   min.segment.length = .1,
                   point.padding=.5,
                   force=5,
                   xlim  = c(x_val_cutoff, NA))+
  
  
  geom_label_repel(data = diff4%>%filter(!is.na(manual_highlights)),
                   aes(label=glabel),
                   color="black",
                   segment.color="black",
                   fill="#faa687",
                   max.overlaps = 100, 
                   point.padding=.5,
                   min.segment.length = .1,
                   force=5)+
                   #xlim  = c(x_val_cutoff, NA))
  
  
  labs(x = expression(~Log[2]~"FC"),
       y = expression("-"~Log[10]~"p-value"),
       title=title_of_plot,
       subtitle = subtitle_of_plot)+
  
  
  theme_minimal()+
  theme(legend.position = "none",
        plot.title = element_text(size = 18, face = "bold"), 
        plot.subtitle = element_text(size=15),# Increase title size
        axis.title = element_text(size = 12))  # Increase axis label size





#vol1


cat("generating grid...\n")

#save plot
plot_object_name_for_session<-"vol1"
plot_filename_output<-paste0(save_location, 
                             "vol-",
                             str_trunc(title_of_plot, 25), "-", 
                             str_trunc(subtitle_of_plot, 12), 
                             output_file_suffix,
                             ".pdf")

metadata_text <- paste0("src: ", 
                        rstudioapi::getSourceEditorContext()$path %>%
                          sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.) %>%
                          sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.), 
                        "\n", "time: ", 
                        lubridate::round_date(Sys.time(), "second"))


text_grob <- grid::textGrob(metadata_text, x=.05, just="left", gp = grid::gpar(fontsize = 9, col = "gray30"))

p3src<-gridExtra::grid.arrange(get(plot_object_name_for_session), text_grob, ncol = 1, heights = c(3, 0.3))



cat("saving plot...\n")



ggsave(filename=plot_filename_output,
       
       plot=p3src,
       limitsize = FALSE,
       
       height=7,
       width=9,
       scale = 1.25,
       dpi=600,
       
       title=paste0(get(plot_object_name_for_session)$labels$title,
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ))




cat("plot saved as ", blue(plot_filename_output), "\n")



cat(green("done\n"))



}
