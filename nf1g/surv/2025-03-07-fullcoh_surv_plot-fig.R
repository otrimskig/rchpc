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
coh1<-readRDS("nf1g/surv/cohorts-2025-03-07-v2.rds")


df1<-coh1%>%
  filter(is.na(exclude))%>%
  filter(!is.na(include_in_surv))



df2<-df1

#read in colors mapping
source("nf1g/colors_map_create.R") #source to get most updated.

#read in colors mapping
col_map<-readRDS("nf1g/surv/colors_map_surv.rds")

#set aspect ratio of plot.
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





curve_category_names<-df1%>%
  select(sym(!!split_curves_by))%>%unique()%>%
  pull()

num_categories <- length(curve_category_names)


plot_subset_values<-df1%>%
  select(!!sym(plot_subset_by))%>%
  unique()%>%
  pull()



split_dfs<-df2 %>%
  group_split(!!sym(split_plots_by))%>%
  
  setNames(sort(unique(df2[[split_plots_by]])))
# split_dfs<-split_dfs[5]

plots<-list()






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
                    "subset: ", plot_subset_clean, ": ", 
                    
                    plot_subset_values),
       
       x = "Days Post Injection",
       y = "% Survival",
       color=split_curves_clean
       
  )+
  theme(plot.title = element_text(hjust = 0),
        aspect.ratio=aspectratio,
        plot.margin = margin(5, 5, 40, 5))






#plots[[1]]



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
  filter(p_value<=.05)%>%
  mutate(group_a = factor(group_a, levels = unique(group_a)),
         group_b = factor(group_b, levels = unique(group_b)))%>%
  mutate(row_id = row_number())%>%


# Create the comp "plot"
  ggplot(.) +
  geom_tile(aes(x = 1, y = row_id, fill = group_a), width = 0.4, height = 0.9, color="#5c5c5c") +
  geom_tile(aes(x = 1.45, y = row_id, fill = group_b), width = 0.4, height = 0.9, color="#5c5c5c") +
  geom_tile(aes(x = 2.5, y = row_id), fill = "white", width = 1, height = 0.9) +  # Empty white tile
  
  # Add p-value text on top of the white tile
  # geom_text(aes(x = 2, y = row_id, label = round(p_value, 5), 
  #               
  #               fontface = ifelse(p_value <= 0.05, "bold", "plain"),
  #               
  #               color = ifelse(p_value <= 0.05, "blue", "black")),
  #           size=2) +
  
  geom_text(aes(x = 2.1, y = row_id, label = format(p_value, scientific = TRUE, digits = 3), #set rounding 
                
                fontface = ifelse(p_value <= .05, "bold", "plain"),
                
            ),
            size=2) +
  
  #scale_color_manual(values = c("blue" = "blue", "black" = "black"))+
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

#comp_plot


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

  


for(i in 1:length(plots)){

  
   
  
  
ggsave(paste0("nf1g/surv/pub/surv_full_cohort_37v2.pdf"),
 
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



