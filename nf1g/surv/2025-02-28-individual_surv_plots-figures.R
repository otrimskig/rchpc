source("libs.R")
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(grid)
library(ggkm)
library(scales)
library(survival)
library(survminer)
library(RColorBrewer)


########dataset read in and construction######################
#reading in current dataset. 
coh1<-readRDS("nf1g/surv/cohorts-2025-02-27.rds")%>%
  mutate(aod = as.numeric(aod))%>%
  mutate(aod=if_else(event==0, 150,aod))



df1<-coh1%>%
  filter(is.na(exclude)|exclude>3)%>%

 filter(is.na(exclude_from_hist_reason)|
          grepl("^fd$", exclude_from_hist_reason)|
          grepl("^fd ", exclude_from_hist_reason)|
          grepl("delay", exclude_from_hist_reason)|
          grepl("^fd$", g_macro_obs)|
          grepl("not collected", exclude_from_hist_reason))%>%
  
  mutate(exclude_hist=if_else(!is.na(exclude_from_hist_reason), "1", NA_character_))%>%
  
  relocate(exclude_hist, exclude, metadata, exclude_from_hist_reason)














te<-coh1%>%
  filter(is.na(exclude)|exclude>3)%>%
  relocate(aod, .after="resultant_geno")%>%
  anti_join(df1, by="mouse_num")%>%
  arrange(resultant_geno, aod)%>%
  
  mutate(exclude=if_else(mouse_num=="25499"|mouse_num=="25500"|mouse_num=="25502", "2", exclude))%>%
  mutate(metadata=if_else(mouse_num=="25499"|mouse_num=="25500"|mouse_num=="25502", 
                         
                           if_else(is.na(metadata), #handle positive results of above in 2 different ways
                                   "non tumor-related reason for death", paste0(metadata, "; non tumor-related reason for death")),
                          
                          metadata)
         )%>%
  
  
  
  mutate(metadata=if_else(grepl("cell line", exclude_from_hist_reason), 
                        
                        if_else(is.na(metadata), #handle positive results of above in 2 different ways
                                "samples used for cell lines", paste0(metadata, "; samples used for cell lines")),
                        
                        metadata))%>%
  
  mutate(exclude=if_else((is.na(exclude)|as.numeric(exclude)>3)&grepl("used for cell lines", metadata), "7", exclude))




  
  # 
  # mutate(metadata=if_else(is.na(metadata), "poor data handling. conflicting data sources or data missing.", metadata))%>%
  # 
  # mutate(exclude=if_else(is.na(exclude)|as.numeric(exclude>3), "3", exclude))



write_csv(te, "nf1g/surv/man_exclusion_check0.csv")



te2<-read_csv("nf1g/surv/man_exclusion_check1 - man_exclusion_check1.csv")%>%
  janitor::clean_names()%>%
  select(mouse_num, exclude, exclude_hist, metadata)%>%
  mutate_at(vars(1:last_col()), ~as.character(.))





coh2<-coh1%>%
  mutate(exclude_hist=NA_character_)%>%
  rows_update(df1, by="mouse_num")%>%
  rows_update(te2, by="mouse_num")%>%
  relocate(mouse_num, exclude_hist, exclude, metadata, exclude_from_hist_reason)
  








ggplot(te, aes(x=resultant_geno, y=aod, color=hist_cat_name))+
  geom_point()+
  theme_minimal()+
  ggrepel::geom_label_repel(aes(label=exclude_from_hist_reason, color=resultant_geno))




col_map<-readRDS("nf1g/surv/colors_map_surv.rds")


aspectratio<-.6






#Plot 01###################################################
#set up splits. grouping variables and splitting variables. 

#set name of actual variable to split on.
split_plots_by<-"full_cohort"

#set how the variable will be printed in the plot. 
split_plots_clean<-"full cohort"


#set how curves will be split.
split_curves_by<-"resultant_geno"

#set how the variable will be printed in the plot. 
split_curves_clean<-"Resultant Geno"



plot_subset_by<-"full_cohort"

plot_subset_clean<-"all"










#get all unique values to assign colors to.
#this will ensure consistent colors throughout these plots.
curve_category_names<-df1%>%
  select(sym(!!split_curves_by))%>%unique()%>%
  pull()

num_categories <- length(curve_category_names)
  
# Choose a color palette based on the number of categories
# colors <- brewer.pal(min(num_categories, 12), "Set1") 

# # Create the named list to assign colors to variables.
# colors_map <- setNames(colors, curve_category_names)
colors_map <- setNames(hue_pal()(num_categories), curve_category_names)

#set up dataset split based on previously decided category.



plot_subset_values<-df1%>%
  select(!!sym(plot_subset_by))%>%
  unique()%>%
  pull()



plot_subset_values<-plot_subset_values[1]



# for (subset in 1:length(plot_subset_values)){
  
  
  df2<-df1%>%
 
     filter(!!sym(plot_subset_by)==plot_subset_values[subset])


          



split_dfs<-df2 %>%
  group_split(!!sym(split_plots_by))%>%
  
  setNames(sort(unique(df2[[split_plots_by]])))
# split_dfs<-split_dfs[5]

plots<-list()





#now run plots based on elements of split. 
for(plot_split_index in 1:length(split_dfs)){
plot_split_index<-1  


spec_df <- as_tibble(split_dfs[[plot_split_index]]) %>%
  setNames(gsub("^d\\$", "", names(as_tibble(split_dfs[[plot_split_index]]))))


plots[[plot_split_index]]<-ggplot(spec_df)+
  geom_km(aes(time = aod, 
              color=!!sym(split_curves_by),
              status = event),
          
          linewidth=2)+
  
  
  facet_wrap(
    vars(!!sym(split_plots_by)),
    labeller = labeller(
      .default = function(x) paste0(split_plots_clean, ": ", x)
    )
  ) +
  
  xlim(0,150)+
  ylim(0,1)+
  
  scale_color_manual(values = col_map[[split_curves_by]])+
  
  theme_classic()+
  #theme(legend.position = c(0.2, 0.2))
  
  labs(title=paste0("Overall Survival ", "\n",
                    "by: ",  split_curves_clean, "\n",
                    "facet: ", split_plots_clean, "\n",
                    "subset: ", plot_subset_clean, ": ", plot_subset_values[subset]),
       x = "Days Post Injection",
       y = "% Survival",
       color=split_curves_clean
       
  )+
  theme(plot.title = element_text(hjust = 0),
        aspect.ratio=aspectratio,
        plot.margin = margin(5, 5, 40, 5))






# plots[[1]]



spec_df_filtered <- spec_df %>%
  group_by(across(all_of(split_curves_by))) %>%
  filter(n() >= 2) %>%
  ungroup()


results <- pairwise_survdiff(
  reformulate(split_curves_by, response = "Surv(time = aod, event = event)"), 
  data = spec_df_filtered, 
  p.adjust.method = "none"
)
    

pvals_table<-results[["p.value"]]
  
if(dim(pvals_table)[1]==0){
    
    print("no comparisons in table.")
    
}else{
    

comps<-results[["p.value"]]%>%
      as_tibble()%>%
      mutate(group_a = attributes(results[["p.value"]])[["dimnames"]][[1]])%>%
      relocate(group_a)%>%
      arrange(group_a)%>%
      pivot_longer(cols=2:last_col(), names_to = "group_b", values_to = "p_value")%>%
      
      mutate(padj_method=results[["p.adjust.method"]])%>%
      
      mutate(grouping_category=results[["grouping"]])%>%
     
      as_tibble()%>%
      filter(!is.na(p_value))
    
   

comp_plot<-comps %>%
  filter(p_value<1)%>%
  mutate(group_a = factor(group_a, levels = unique(group_a)),
         group_b = factor(group_b, levels = unique(group_b)))%>%
  mutate(row_id = row_number())%>%


# Create the plot
  ggplot(.) +
  geom_tile(aes(x = 1, y = row_id, fill = group_a), width = 0.4, height = 0.9, color="#5c5c5c") +
  geom_tile(aes(x = 1.45, y = row_id, fill = group_b), width = 0.4, height = 0.9, color="#5c5c5c") +
  geom_tile(aes(x = 2.5, y = row_id), fill = "white", width = 0.4, height = 0.9) +  # Empty white tile
  
  # Add p-value text on top of the white tile
  # geom_text(aes(x = 2, y = row_id, label = round(p_value, 5), 
  #               
  #               fontface = ifelse(p_value <= 0.05, "bold", "plain"),
  #               
  #               color = ifelse(p_value <= 0.05, "blue", "black")),
  #           size=2) +
  
  geom_text(aes(x = 2, y = row_id, label = format(p_value, scientific = TRUE, digits = 7), 
                
                fontface = ifelse(p_value <= 10, "bold", "plain"),
                
            ),
            size=4) +
  
  scale_color_manual(values = c("blue" = "blue", "black" = "black"))+
  scale_fill_manual(values = col_map[[split_curves_by]]) +
  scale_x_continuous(breaks = c(1, 1.45, 2), labels = c("Group A", "Group B", "p-value")) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none")+
  # 
  coord_fixed(ratio = .2)

comp_plot


# combined_plot <- grid.arrange(
#   grobs = list(plots[[1]], comp_plot),
#   nrow = 1,  # Arrange the plots in a single column
#   widths = c(10, 1)  # Set the heights relative to each plot
# )

plots[[plot_split_index]]<-plots[[plot_split_index]] +
  annotation_custom(
    grob = ggplotGrob(comp_plot),  # The grob of the second plot
    xmin = 0, xmax = 40,  # Set the x-range for placing the plot inside
    ymin = 0, ymax = .5   # Set the y-range for placing the plot inside
  )

}
}














for(i in 1:length(plots)){
  
  metadata_text <- paste0("src: ", 
                          rstudioapi::getSourceEditorContext()$path %>%
                            sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.) %>%
                            sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.), 
                          "\n", "time: ", 
                          lubridate::round_date(Sys.time(), "second"))
  
  text_grob <- textGrob(metadata_text, x=.05, just="left", gp = gpar(fontsize = 7, col = "gray30"))
  
 plots[[i]]<-grid.arrange(plots[[i]], text_grob, ncol = 1, heights = c(3, 0.3))
  
  
  
   
   
   
}   

  
  
  
  
stop()
for(i in 1:length(plots)){

  
   
  
  
ggsave(paste0("nf1g/surv/plots/surv_", split_plots_by, "-", split_curves_by, "-", plot_subset_clean, "-plot_", subset, "-", i,".pdf"),
 
        title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       plot=plots[[i]],
       limitsize = FALSE,
       
       
       height=5,
       width=8,
       scale = 1.2,
       dpi=600,
       
       
       
)
}


}




















































































plot_paths<-fs::dir_info("nf1g/surv/plots", full.names = T, pattern="*\\.pdf$", recurse = F)%>%tibble()%>%
  filter(modification_time>Sys.time()-lubridate::minutes(1))%>%
  arrange(path)%>%
  filter(type=="file")%>%
  pull(path)




combined_output_filen<-paste0("nf1g/surv/plots/comb/surv_comb_", split_plots_by, "-", split_curves_by, "-", plot_subset_clean, ".pdf")



qpdf::pdf_combine(plot_paths, output = combined_output_filen)




