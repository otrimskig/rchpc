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

#read in colors mapping
source("nf1g/colors_map_create.R") #source to get most updated.
col_map<-readRDS("nf1g/surv/colors_map_surv.rds")

plot_output_loc<-"nf1g/surv/pub/pub_plots/surv_all/" #include "/"
plot_output_base<-"surv_looped_" #do not include .pdf


########dataset read in and construction############################
#reading in current dataset. 
coh1<-readRDS("nf1g/surv/cohorts-2025-03-07-v2.rds")

df1<-coh1%>%
  filter(is.na(exclude))%>%
  filter(!is.na(include_in_surv))
# filter(is.na(exclude_hist))


#set aspect ratio of plot.
aspectratio<-.6

#Plot 01############################################################
#set up splits. grouping variables and splitting variables. 

#set name of actual variable to split on.
split_plots_by<-"hist_cat_name"

#set how the variable will be printed in the plot. 
split_plots_clean<-"Tumor Type (Histology)"

#set how curves will be split.
split_curves_by<-"resultant_geno"

#set how the variable will be printed in the plot. 
split_curves_clean<-"Cohort (Resultant Genotype)"

#
# plot_subset_by<-"hist_grade_name"
# # 
# plot_subset_clean<-"Histology Grade"

plot_subset_by<-"hist_cat_name"
# 
plot_subset_clean<-"Tumor Type (Histology)"

su

#anything else that can be added to df. #############################
plot_subset_values<-df1%>%
  select(!!sym(plot_subset_by))%>%
  unique()%>%
  pull()%>%
  as.character()


# plot_subset_values<-plot_subset_values[1]
for (su in 1:length(plot_subset_values)){
  print(paste0("working on: su: ",su)) 
  subset<-as.character(plot_subset_values[su])

if(!is.na(subset)){
  
  df2<-df1%>%
    filter(!!sym(plot_subset_by)==subset)%>%
    mutate(across(where(is.factor), droplevels))

      
print(df2)
print(subset)

if(nrow(df2)>0){

#Set up dfs and initialize plot list#################################
curve_category_names<-df2%>%
  select(sym(!!split_curves_by))%>%unique()%>%
  pull()

num_categories <- length(curve_category_names)

plot_subset_values<-df2%>%
  select(!!sym(plot_subset_by))%>%
  unique()%>%
  pull()

split_dfs<-df2 %>%
  group_split(!!sym(split_plots_by))%>%
  setNames(sort(unique(df2[[split_plots_by]])))
# split_dfs<-split_dfs[5]

plot_list_v1<-list()


#create initial plots################################################
for (plot_split_index in 1:length(split_dfs)){

  spec_df1 <- as_tibble(split_dfs[[plot_split_index]]) %>%
    setNames(gsub("^d\\$", "", names(as_tibble(split_dfs[[plot_split_index]]))))
  
  cat_name<-spec_df1%>%
    select(all_of(split_plots_by))%>%
    mutate_at(vars(all_of(split_plots_by)), ~as.character(.x))%>%
    unique()%>%
    pull()
  
  spec_counts0<-reframe(spec_df1, .by=all_of(split_curves_by),
                        countn=n())
  
  spec_counts1 <- spec_df1 %>%
    select(all_of(split_curves_by)) %>%      # Select column by name dynamically
    distinct()%>%# Get unique rows
    left_join(spec_counts0, by = split_curves_by) %>%  # Dynamically use the split_curves_by for joining
    arrange(!!sym(split_curves_by)) %>%      # Dynamically use the split_curves_by for arranging
    mutate(var_plus_countn = paste0(         # Concatenate values for var_plus_countn
      .[[split_curves_by]], " (n = ", countn, ")"          # Reference the column dynamically using `[[` for concatenation
    )) %>%
    mutate(var_plus_countn = factor(var_plus_countn, 
                                    levels = var_plus_countn[order(!!sym(split_curves_by))]))
  
  
  
  plot_list_v1[[plot_split_index]]<-ggplot(spec_df1)+
    geom_km(aes(time = aod, 
                color=!!sym(split_curves_by),
                status = event),
    linewidth=2.5
    # alpha= don't set here. geom_km doesn't respect. set in scale_color_manual.
    )+
    
    facet_wrap(
      vars(!!sym(split_plots_by)),
      labeller = labeller(
        .default = function(x) paste0(split_plots_clean, ": ", x)
      )
    ) +
    
    xlim(0,150)+
    ylim(0,1)+
    
    scale_color_manual(values = alpha(col_map[[split_curves_by]], 0.9),  #geom_km doesn't respect alpha argument for some reason. set it manually in color call.
                       labels = spec_counts1$var_plus_countn)+ 
    theme_classic()+
    
    labs(title=paste0("Overall Survival ", "\n",
                      "by: ",  split_curves_clean, "\n",
                      "facet: ", split_plots_clean, "\n",
                      "subset: ", plot_subset_clean, ": ", 
                      
                      cat_name),
         
         x = "Days Post Injection",
         y = "% Survival",
         color=split_curves_clean
         
    )+
    theme(plot.title = element_text(hjust = 0),
          aspect.ratio=aspectratio,
          plot.margin = margin(5, 5, 5, 5))

} #end of initial plot create.






plot_list_v1[[1]]


plot_list_v2<-list()














#create grob for thing###############################################
# for (plot_split_index in 1:length(split_dfs)){
#   
#   spec_df1 <- as_tibble(split_dfs[[plot_split_index]]) %>%
#     setNames(gsub("^d\\$", "", names(as_tibble(split_dfs[[plot_split_index]]))))
#   
#   spec_df_filtered <- spec_df1 %>%
#     group_by(across(all_of(split_curves_by))) %>%
#     filter(n() >= 2) %>%
#     ungroup()
#   
#   results <- pairwise_survdiff(
#     reformulate(split_curves_by, response = "Surv(time = aod, event = event)"), 
#     data = spec_df_filtered, 
#     p.adjust.method = "none"
#   )
#   
#   pvals_table<-results[["p.value"]]
#   
#   if(dim(pvals_table)[1]==0){
#     
#     print("no comparisons in table.")
#     
#   }else{

for (plot_split_index in 1:length(split_dfs)) {
  
  spec_df1 <- as_tibble(split_dfs[[plot_split_index]]) %>%
    setNames(gsub("^d\\$", "", names(as_tibble(split_dfs[[plot_split_index]]))))
  
  spec_df_filtered <- spec_df1 %>%
    group_by(across(all_of(split_curves_by))) %>%
    filter(n() >= 2, any(event == 1)) %>%  # Ensure at least one event in each group
    ungroup()
  
  if (n_distinct(spec_df_filtered[[split_curves_by]]) >= 2) {
    
  results <- pairwise_survdiff(
    reformulate(split_curves_by, response = "Surv(time = aod, event = event)"), 
    data = spec_df_filtered, 
    p.adjust.method = "none"
  )
  
  pvals_table <- results[["p.value"]]
  
  if (dim(pvals_table)[1] > 0) {

    comps<-results[["p.value"]]%>%
      as_tibble()%>%
      mutate(group_a = attributes(results[["p.value"]])[["dimnames"]][[1]])%>%
      relocate(group_a)%>%
      arrange(group_a)%>%
      pivot_longer(cols=2:last_col(), names_to = "group_b", values_to = "p_value")%>%
      
      mutate(padj_method=results[["p.adjust.method"]])%>%
      
      mutate(grouping_category=results[["grouping"]])%>%
      
      as_tibble()%>%
      filter(!is.na(p_value))%>%
      filter(p_value<=.05)
    
    if(dim(comps)[1]==0){ #check for dimensions in comps tibble (no significant comps)
      
      #make a no comps are significant square annotation. 
      comp_plot<-ggplot(comps) +
        
        geom_tile(aes(x = 7.5, y = 7.5), fill = "#ccffff", alpha=0, width = 15, height = 15)+  # Empty tile
        geom_tile(aes(x = .75, y = 1, fill = "gray"), width = 1.5, height = .85, color="#5c5c5c")+
        geom_tile(aes(x = 2.35, y = 1, fill = "gray"), width = 1.5, height = .85, color="#5c5c5c")+
        
        
        geom_text(aes(x = 3.3, y = 1, label = "all comparisons ns",
                      
                      fontface ="bold"
        ),
        size=2,
        hjust = 0) +  
        
        #scale_fill_manual(values = col_map[[split_curves_by]]) +
        theme_minimal() +
        
        theme(
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none",
          plot.margin = margin(0, 0, 0, 0))+
        coord_fixed(ratio = 1)
      
    }else{#gray plot 
    }else{
      
      
      comp_plot<-comps %>%
        filter(p_value<=.05)%>%
        mutate(group_a = factor(group_a, levels = unique(group_a)),
               group_b = factor(group_b, levels = unique(group_b)))%>%
        mutate(row_id = row_number())%>%
        
        
        # Create the comp "plot"
        ggplot(.) +
        geom_tile(aes(x = 7.5, y = 7.5), fill = "#ccffff", alpha=0, width = 15, height = 15)+  # Empty tile
        geom_tile(aes(x = .75, y = row_id, fill = group_a), width = 1.5, height = .85, color="#5c5c5c")+
        geom_tile(aes(x = 2.35, y = row_id, fill = group_b), width = 1.5, height = .85, color="#5c5c5c")+
        
        
        geom_text(aes(x = 3.3, y = row_id, label = format(p_value, scientific = TRUE, digits = 3), #set rounding
                      
                      fontface = ifelse(p_value <= .05, "bold", "plain"),
                      
        ),
        size=2,
        hjust = 0) +
        
        scale_fill_manual(values = col_map[[split_curves_by]]) +
        theme_minimal() +
        
        theme(
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none",
          plot.margin = margin(0, 0, 0, 0))+
        coord_fixed(ratio = 1)
      
      
    }
  }
    
    plot_list_v2[[plot_split_index]]<-plot_list_v1[[plot_split_index]] +
      annotation_custom(
        grob = ggplotGrob(comp_plot),  # The grob of the second plot
        xmin = -10, xmax = 40,  # Set the x-range for placing the plot inside
        ymin = 0, ymax = .5   # Set the y-range for placing the plot inside
      )
    
    
  }
}
# 
# if(length(plot_list_v2)>=0){}

plot_list_v3<-list()
#add in meta src data below plots.####################################
for(i in 1:length(plot_list_v2)){
  
  metadata_text <- paste0("src: ", 
                          rstudioapi::getSourceEditorContext()$path %>%
                            sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.) %>%
                            sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.), 
                          "\n", "time: ", 
                          lubridate::round_date(Sys.time(), "second"))
  
  text_grob <- textGrob(metadata_text, x=.05, just="left", gp = gpar(fontsize = 7, col = "gray30"))
  
  plot_list_v3[[i]] <- grid.arrange(
    plot_list_v2[[i]], 
    text_grob, 
    ncol = 1, 
    heights = c(3, 0.3),  # Maintain spacing
    layout_matrix = rbind(c(1), c(2))  # Keeps layout structure stable
  )

}


#grid.draw(plot_list_v3[[1]])

##save all plots to pdfs############################################
for(i in 1:length(plot_list_v3)){

  ggsave(paste0(plot_output_loc,
                plot_output_base,
                split_plots_by, "-", split_curves_by, "-", plot_subset_clean, "-plot_", subset, "-", i,".pdf"),

         title=paste0(gsub("\n", "][", plot_list_v1[[1]]$labels$title), "][",

           "src: ",
                      rstudioapi::getSourceEditorContext()$path%>%
                        sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                        sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                      
                      " at ", 
                      
                      lubridate::round_date(Sys.time(), "second")
         ),
         
         plot=plot_list_v3[[i]],
         limitsize = FALSE,
         
         height=5,
         width=10,
         scale = 1.2,
         dpi=600,
  )
}

}else{
  print("no data in table. going to next")
}

}else{
  print("na subset.")
}

}

