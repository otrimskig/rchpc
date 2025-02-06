source("libs.R")
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(ggkm)
library(scales)
library(survival)
library(survminer)
library(RColorBrewer)


########dataset read in and construction######################
#reading in current dataset. 
coh1<-readRDS("nf1g/surv/cohorts-2025-01-07.rds")%>%
  mutate(aod = as.numeric(aod))

#final form for survival analysis. 
df0<-coh1%>%
  
  filter(is.na(exclude_from_hist_reason))%>%
  
  
  #keep only necessary vars
  select(exclude, strain_injection, resultant_geno, 
         
          hist,
         
         
         aod, exp_end_date)%>%
  
  mutate(event = if_else(is.na(exp_end_date),1,0))%>%
  select(-exp_end_date)%>%
  
  #making a convenient variable to easily compare nf1 vs. non-nf1s. 
  mutate(nf1=substr(resultant_geno, 1, 6))%>%
  mutate(geno_bg=substr(resultant_geno, 8, nchar(resultant_geno)))%>%
  
  #most conservative exclusion criteria
  filter(is.na(exclude)|exclude>3)%>%
  
  mutate(aod=if_else(event==0, 150,aod))%>%
  
  mutate(hist_cat=substr(hist, 1,1))%>%
  mutate(hist_grade=substr(hist, 3,4))%>%
  mutate(hist_grade=gsub("\\.", "", hist_grade))%>%
  mutate(hist_grade=gsub("15", "1.5", hist_grade))%>%
  mutate(hist_grade=gsub("^d$", "ned", hist_grade))%>%
  
  
  
  mutate(hist_catgrade=paste0(hist_cat, ".", hist_grade))%>%
  mutate(hist_catgrade=if_else(hist_catgrade=="n.ned", "ned", hist_catgrade))%>%
  
  mutate(hist_cat=if_else(hist_cat=="n", "ned", hist_cat))%>%
  mutate(hist_grade=if_else(hist_grade=="N", "no grade assigned", hist_grade))



nf1_renaming_dataset <- readxl::read_excel("nf1g/surv/nf1-renaming dataset.xlsx", 
                                   sheet = "patho_cat")

df1<-df0%>%
  left_join(nf1_renaming_dataset)


aspectratio<-.6


#make empty list
plot_list <- list() 
master_plot_list <- list() 
#overall survival: all cohorts. 




#Plot 01###################################################
#set up splits. grouping variables and splitting variables. 

#set name of actual variable to split on.
split_plots_by<-"hist_grade"

#set how the variable will be printed in the plot. 
split_plots_clean<-"Histology Grade"


#set how curves will be split.
split_curves_by<-"patho_cat_name"

#set how the variable will be printed in the plot. 
split_curves_clean<-"Tumor Category"

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

split_dfs<-df1 %>%
  group_split(!!sym(split_by))%>%
  
  setNames(sort(unique(df1[[split_by]])))
split_dfs<-split_dfs[5]

plots<-list()

#now run plots based on elements of split. 
for(plot_split_index in 1:length(split_dfs)){

spec_df <- as_tibble(split_dfs[[plot_split_index]]) %>%
  setNames(gsub("^d\\$", "", names(as_tibble(split_dfs[[plot_split_index]]))))


plots[[plot_split_index]]<-ggplot(spec_df)+
  geom_km(aes(time = aod, 
              color=patho_cat_name,
              status = event),
          
          linewidth=2)+
  
  
  facet_wrap(vars(hist_grade), labeller = labeller(hist_grade = function(x) paste(split_plots_clean, x)))+
  # 
  
  xlim(0,150)+
  
  scale_color_manual(values = colors_map)+
  
  theme_classic()+
  #theme(legend.position = c(0.2, 0.2))
  
  labs(title=paste0("Overall Survival ", "\n",
                    "by: ",  split_curves_clean, "\n",
                    "facet: ", split_plots_clean, "\n",
                    "subset: ", plot_subset_clean),
       x = "Days Post Injection",
       y = "% Survival",
       color=split_curves_clean
       
  )+
  theme(plot.title = element_text(hjust = 0),
        aspect.ratio=aspectratio)

}



plots[1]



spec_df


results <- pairwise_survdiff(Surv(time = aod, event = event) ~ patho_cat_name, 
                                data = spec_df, 
                                p.adjust.method = "none")
    

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
    
}   



comps


comp_plot<-comps %>%
  mutate(group_a = factor(group_a, levels = unique(group_a)),
         group_b = factor(group_b, levels = unique(group_b)))%>%
  mutate(row_id = row_number())%>%


# Create the plot
  ggplot(.) +
  geom_tile(aes(x = 1, y = row_id, fill = group_a), width = 0.4, height = 0.9) +
  geom_tile(aes(x = 1.45, y = row_id, fill = group_b), width = 0.4, height = 0.9) +
  geom_tile(aes(x = 1.9, y = row_id), fill = "white", width = 0.4, height = 0.9) +  # Empty white tile
  
  # Add p-value text on top of the white tile
  geom_text(aes(x = 1.9, y = row_id, label = round(p_value, 3), 
                fontface = ifelse(p_value <= 0.05, "bold", "plain"),
                color = ifelse(p_value <= 0.05, "blue", "black"))) +
  
  
  
  scale_color_manual(values = c("blue" = "blue", "black" = "black"))+
  scale_fill_manual(values = colors_map) +
  scale_x_continuous(breaks = c(1, 1.45, 1.9), labels = c("Group A", "Group B", "p-value")) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )+
  coord_fixed(ratio = .2)




comp_plot

plots[1]



combined_plot <- grid.arrange(
  grobs = list(plots[[1]], comp_plot),
  ncol = 1,  # Arrange the plots in a single column
  heights = c(10, 1)  # Set the heights relative to each plot
)

plots[[1]] +
  annotation_custom(
    grob = ggplotGrob(comp_plot),  # The grob of the second plot
    xmin = 0, xmax = 40,  # Set the x-range for placing the plot inside
    ymin = 0, ymax = .25   # Set the y-range for placing the plot inside
  )
