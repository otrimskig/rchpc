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
library(crayon)
library(ggprism)
library(ggh4x)



#read in colors mapping
source("nf1g/colors_map_create.R") #source to get most updated.
col_map<-readRDS("nf1g/surv/colors_map_surv.rds")

plot_output_loc<-"nf1g/surv/pub/pub_plots/surv_prism/" #include "/"
plot_output_base<-"surv_" #do not include .pdf


########dataset read in and construction############################
#reading in current dataset. 
coh1<-readRDS("nf1g/surv/cohorts-2025-03-07-v2.rds")

df1<-coh1%>%
  filter(is.na(exclude))%>%
  filter(!is.na(include_in_surv))%>%
  
  
  
  slice(1:20)
# filter(is.na(exclude_hist))


#set aspect ratio of plot.
aspectratio<-.6

#Plot 01############################################################
#set up splits. grouping variables and splitting variables. 

#set name of actual variable to split on.
split_plots_1deg<-"hist_cat_name"

#set how the variable will be printed in the plot. 
split_plots_1deg_clean<-"Tumor Type (Histology)"



#set further subsets if desired.
split_plots_2deg<-"hist_grade_name"
# 
split_plots_2deg_clean<-"Histology Grade"



# 
# 
# #set further subsets if desired.
# split_plots_2deg<-"hist_cat_name"
# # 
# split_plots_2deg_clean<-"Tumor Type (Histology)"







#set how curves will be split.
split_curves_by<-"resultant_geno"

#set how the variable will be printed in the plot. 
split_curves_clean<-"Cohort (Resultant Genotype)"

#####gray plot generation################################

#make a no comps are significant square annotation. 
gray_plot<-ggplot(sample_n(mtcars, 1)) +
  
  geom_tile(aes(x = 7.5, y = 7.5), fill = "#ccffff", alpha=0, width = 15, height = 15)+  # Empty tile
  geom_tile(aes(x = .75, y = 1), fill = "gray", width = 1.5, height = .85, color="#5c5c5c")+
  geom_tile(aes(x = 2.35, y = 1), fill = "gray", width = 1.5, height = .85, color="#5c5c5c")+
  
  geom_text(aes(x = 3.3, y = 1), 
            label = "comparison is NA",
            fontface ="bold",
            size=2,
            hjust = 0) +  
  
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


######################################################################



plot_list_v1<-list()
plot_list_v2<-list()


#anything else that can be added to df. #############################
split_1deg_values<-df1%>%
  select(!!sym(split_plots_1deg))%>%
  unique()%>%
  pull()%>%
  as.character()



# split_1deg_values<-split_1deg_values[1]
for (sp1 in 1:length(split_1deg_values)){
  
  
  #working on each split_1_char.
  cat("\n", "\n",bold(blue("Working on: ")), split_1deg_values[[sp1]], "\n")
  
  split_1_char<-as.character(split_1deg_values[sp1])
  
  if(is.na(split_1_char)){
    
    #split_1_char is na. exit.
    cat(bold(red("Is entirely NA:")), split_1_char, "\n")
    
  }else{
    
    #split_1_char exists. make df2.
    cat(bold(cyan("split_1_char exists: ")), split_1_char, "\n")
    
    df2<-df1%>%
      filter(!!sym(split_plots_1deg)==split_1_char)%>%
      mutate(across(where(is.factor), droplevels))
    
    #check df2 contains data.
    if(nrow(df2)==0){  
      
      #no data in df2.
      cat(bold(red("No values in df2: ")), split_1_char, "\n")
      
    }else{
      
      plot_list_v1[[split_1_char]] <- list()
      plot_list_v2[[split_1_char]] <- list()
      
      #data in df2.
      cat(bold(cyan("df2 has values: ")), split_1_char, "\n")
      
      
      #Set up dfs and initialize plot list#################################
      split_dfs<-df2 %>%
        group_split(!!sym(split_plots_2deg))%>%
        setNames(sort(unique(df2[[split_plots_2deg]])))
      
      #loop through split_dfs.
      for (sp2 in 1:length(split_dfs)) {
        
        split_2_char<-names(split_dfs)[sp2]
        
        spec_df1 <- as_tibble(split_dfs[[sp2]]) %>%
          setNames(gsub("^d\\$", "", names(as_tibble(split_dfs[[sp2]]))))
        
        if(nrow(spec_df1)==0){  
          
          #no data in spec_df1.
          cat(bold(red("No values in spec_df1: ")), split_2_char, "\n")
          
        }else{
          
          #sufficient data in specific dataframe to proceed.
          cat(bold(cyan("Sufficient values in spec_df1")), split_2_char, "\n")
          
          #initialize nested lists.
          plot_list_v1[[split_1_char]][[split_2_char]] <- list()
          plot_list_v2[[split_1_char]][[split_2_char]] <- list()
          
          #get counts of each group within specific df.
          spec_counts0<-reframe(spec_df1, .by=all_of(split_curves_by),
                                countn=n())
          
          #join count number to main spec df. 
          #This ensures that the factor order is retained.
          #create the string which includes the group name plus group count for each, that
          #will be used on the plot.
          spec_counts1 <- spec_df1 %>%
            select(all_of(split_curves_by)) %>%      # Select column by name dynamically
            distinct()%>%# Get unique rows
            left_join(spec_counts0, by = split_curves_by) %>%  # Dynamically use the split_curves_by for joining
            arrange(!!sym(split_curves_by)) %>%      # Dynamically use the split_curves_by for arranging
            mutate(var_plus_countn = paste0(         # Concatenate values for var_plus_countn
              .[[split_curves_by]], " (n = ", countn, ")"          # Reference the column dynamically using `[[` for concatenation
            )) %>%
            mutate(var_plus_countn = factor(var_plus_countn, 
                                            levels = var_plus_countn[order(!!sym(split_curves_by))]))#ensure order of new string matches order of factored groups.
          
          # Construct the title outside of labs()
          title_text <- if (split_1_char == split_2_char) {
            paste0("Overall Survival ", "\n",
                   "by: ", split_curves_clean, "\n",
                   "split 1deg: ", split_plots_1deg_clean, ": ", split_1_char, "\n",
                   "split 2deg: None")
          } else {
            paste0("Overall Survival ", "\n",
                   "by: ", split_curves_clean, "\n",
                   "split 1deg: ", split_plots_1deg_clean, ": ", split_1_char, "\n",
                   "split 2deg: ", split_plots_2deg_clean, ": ", split_2_char)
          }
          
          
          
          
          #create plot and assign to list location.
          plot_list_v1[[split_1_char]][[split_2_char]]<-ggplot(spec_df1)+
            geom_km(aes(time = aod, 
                        color=!!sym(split_curves_by),
                        status = event),
                    linewidth=1.5
                    # alpha= don't set here. geom_km doesn't respect. set in scale_color_manual.
            )+
            
            
            # facet_wrap(
            #   vars(!!sym(split_plots_1deg)),
            #   labeller = labeller(
            #     .default = function(x) {
            #       if (split_1_char == split_2_char) {
            #         "poop" #split_1_char
            #       } else {
            #         "pee" #paste0(split_1_char, ": ", split_2_char)
            #       }
            #     }
            #   )
            # ) +  
            # 
            
            scale_x_continuous(
              limits = c(0, 150),
              expand = c(0, 0)
            )+
            
            
            scale_color_manual(values = alpha(col_map[[split_curves_by]], 0.9),  #geom_km doesn't respect alpha argument for some reason. set it manually in color call.
                               labels = spec_counts1$var_plus_countn)+ 
            theme_classic() +
            
            # Pass the constructed title to labs()
            labs(
              title = title_text,
              x = "Days Elapsed",
              y = "Probability of Survival",
              color = split_curves_clean
            ) +
            
            
            
            
            
            theme_prism() +
            
            guides(
              y = guide_axis(minor.ticks =
                               TRUE),
              color = guide_legend(ncol = 1, byrow = TRUE) # Force a line break after each item) +
            )+
            
            
            scale_y_continuous(
              limits = c(0, 1.1),
              breaks = c(0, 0.5, 1),
              labels = c("0", "50", "100"),
              minor_breaks = seq(0, 1.1, by = 0.05),
              expand = c(0, 0)
            )+
            
            
            theme(,
                  aspect.ratio=aspectratio,
                  plot.margin = margin(20, 20, 20, 20),
                  axis.ticks.length = unit(10, "pt"),         # Make major ticks longer
                  axis.minor.ticks.length = unit(5, "pt"),
                  axis.text.y = element_text(
                    margin = margin(r = 4)),
                  plot.title = element_text(hjust = 0, margin = margin(b = 50)),
                  
                  
                  legend.margin = margin(t = 20),
                  legend.position = "bottom",
                  legend.justification = "left",
                  
                  #adjust size of legend txt
                  legend.text = element_text(size = 12)
            )
          
          ##############################
          #end of initial plot creation
          ############################
          
          
          
          
          #now run tests to and either make non-sig or significant annotation.
          if (n_distinct(spec_df1[[split_curves_by]]) <2) {
            
            #bad distinction between values.
            cat(bold(red("not enough distinct values in: ")), split_1_char, ":", split_2_char, "\n")
            
            #set comp_plot to gray plot.
            comp_plot<-gray_plot
            
            
          }else{
            
            #enough distinction to proceed with p-value assessment.
            cat(bold(cyan("enough distinct values in: ")), split_2_char, "\n")
            
            
            results <- pairwise_survdiff(
              reformulate(split_curves_by, response = "Surv(time = aod, event = event)"), 
              data = spec_df1, 
              p.adjust.method = "none"
            )
            
            pvals_table <- results[["p.value"]]
            
            
            if (dim(pvals_table)[1] == 0) { 
              
              #bad distinction in p value table
              cat(bold(red("not enough distinct values in p value table: ")), split_1deg_values[[sp1]], ":", split_1_char, 
                  ":", sp2,
                  
                  "\n")
              
              #set comp_plot to gray plot.
              comp_plot<-gray_plot
              
            }else{
              
              #enough distinction in p value table.
              cat(bold(cyan("sufficient values in p value table: ")), split_2_char, "\n")
              
              comps<-results[["p.value"]]%>%
                as_tibble()%>%
                mutate(group_a = attributes(results[["p.value"]])[["dimnames"]][[1]])%>%
                relocate(group_a)%>%
                arrange(group_a)%>%
                pivot_longer(cols=2:last_col(), names_to = "group_b", values_to = "p_value")%>%
                
                mutate(padj_method=results[["p.adjust.method"]])%>%
                
                mutate(grouping_category=results[["grouping"]])%>%
                
                as_tibble()%>%
                filter(!is.na(p_value))#%>%
              #filter(p_value<=.05)
              
              if(dim(comps)[1]==0){
                #no p values generated
                cat(bold(red("no p values generated: ")), split_1_char, ":", split_2_char, "\n")
                
                #set comp_plot to gray plot.
                comp_plot<-gray_plot
                
                
              }else{
                
                #enough distinction in p value table.
                cat(bold(green("results in p value table: ")), split_1_char, ": ", split_2_char, "\n")
                
                comp_plot <- comps %>%
                  arrange(p_value) %>%
                  mutate(
                    group_a = factor(group_a, levels = unique(group_a[order(p_value)])),
                    group_b = factor(group_b, levels = unique(group_b[order(p_value)])),
                    row_id = row_number()
                  )%>%
                  
                  ggplot() +
                  geom_tile(aes(x = 7.5, y = 7.5), fill = "#ccffff", alpha = 0, width = 15, height = 15) +  # Empty tile
                  geom_tile(aes(x = 0.75, y = row_id, fill = group_a), width = 1.5, height = 0.85, color = "#5c5c5c") +
                  geom_tile(aes(x = 2.35, y = row_id, fill = group_b), width = 1.5, height = 0.85, color = "#5c5c5c") +
                  
                  geom_text(aes(
                    x = 3.3, y = row_id, 
                    label = format(p_value, scientific = TRUE, digits = 3),
                    fontface = ifelse(p_value <= 0.05, "bold", "plain"),
                    color = ifelse(p_value <= 0.05, "black", "gray50")
                  ), size = 2, hjust = 0) +
                  
                  scale_color_identity()+
                  
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
          }}}}}}
          
          
          #after this, comp_plot should be defined as either a gray
          #or a generated sig plot.
          
          #add in p value annotation.
plot_list_v1[[split_1_char]][[split_2_char]]+
            
            annotation_custom(
              grob = ggplotGrob(comp_plot)  # The grob of the second plo
            )



plot_list_v1[[split_1_char]][[split_2_char]]+
  
  comp_plot










