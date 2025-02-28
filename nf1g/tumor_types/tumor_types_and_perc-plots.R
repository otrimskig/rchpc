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
coh1<-readRDS("nf1g/surv/cohorts-2025-02-27.rds")

df1<-coh1%>%
#most conservative exclusion criteria
filter(is.na(exclude)|exclude>3)


aspectratio<-.6


genotypes<-df1%>%select(resultant_geno)%>%unique()%>%arrange()%>%pull()





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
  left_join(coh1%>%select(resultant_geno, resultant_geno_name)%>%unique())



geno_subset$resultant_geno_numeric <- as.numeric(as.factor(geno_subset$resultant_geno))
geno_levels <- unique(geno_subset$resultant_geno)
geno_subset$resultant_geno_numeric <- match(geno_subset$resultant_geno, geno_levels) * 0.6  # Adjust 0.8 to fine-tune spacing

p <- ggplot(geno_subset) +
  geom_col(aes(x = resultant_geno_numeric, 
               fill = patho_cat_name,
               y = perc),
           width = 0.4,  
           position = position_dodge(0.5)) +  
  theme_classic() +

  theme(
    axis.text.x = element_text(size=12,angle = 45, hjust = 1),
    plot.margin = margin(100, 100, 100, 100)
  ) +
  ggtitle("Tumor type prevalence by resultant genotype")+
  labs(fill="Tumor Type",
       x=NULL,
       y="% of Cohort")+
  scale_x_continuous(breaks = unique(geno_subset$resultant_geno_numeric),
                     labels = unique(geno_subset$resultant_geno))  # Keep original labels

p
















  
  
ggsave("nf1g/tumor_types/tumor_types.pdf",
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       plot=p,
       limitsize = FALSE,
       
       
       height=10,
       width=10,
       scale = 1,
       dpi=600,
       
       
       
)  
  






###############################################################





genotypes<-df1%>%select(resultant_geno)%>%unique()%>%arrange()%>%pull()


geno_subset<-df1%>%
  #filter(resultant_geno==genotypes[1]|resultant_geno==genotypes[2])%>%
  select(resultant_geno, patho_cat_name)%>%
  group_by(resultant_geno,patho_cat_name)%>%
  summarise(n=n(),
  )%>%
  ungroup()%>%
  group_by(resultant_geno)%>%
  summarize(patho_cat_name, n=n,total_n=sum(n))%>%
  mutate(perc=n/total_n*100)%>%
  ungroup()%>%
  complete(resultant_geno, patho_cat_name, fill = list(perc = -1))



geno_subset


prop_test_results <- geno_subset %>%
  group_by(patho_cat_name) %>%
  summarise(
    p_value = prop.test(
      x = n,  # The count of successes (e.g., number of "green" individuals)
      n = total_n  # The total sample size
    )$p.value  # Extract p.value directly from the result
  )





color_map<-readRDS("nf1g/ds/colors_map_geno.rds")


geno_subset$patho_cat_name_numeric <- as.numeric(as.factor(geno_subset$patho_cat_name))
cat_levels <- unique(geno_subset$patho_cat_name)
geno_subset$patho_cat_name_numeric <- match(geno_subset$patho_cat_name, cat_levels) * 0.6  # Adjust 0.8 to fine-tune spacing

p2 <- ggplot(geno_subset) +
  geom_col(aes(x = patho_cat_name_numeric, 
               fill = resultant_geno,
               y = perc),
           width = 0.4,  
           position = position_dodge(0.5)) +  
  theme_classic()+
  scale_fill_manual(values=color_map$resultant_geno)+
  
  theme(
    axis.text.x = element_text(size=12,angle = 45, hjust = 1),
    plot.margin = margin(100, 100, 100, 100)
  ) +
  ggtitle("Genotype prevalence by tumor type")+
  labs(fill="Resultant genotype",
       x=NULL,
       y="% of Cohort")+
  scale_x_continuous(breaks = unique(geno_subset$patho_cat_name_numeric),
                     labels = unique(geno_subset$patho_cat_name))  # Keep original labels

p2



color_map$resultant_geno




ggsave("nf1g/tumor_types/tumor_types-2.pdf",
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       plot=p2,
       limitsize = FALSE,
       
       
       height=10,
       width=10,
       scale = 1,
       dpi=600,
       
       
       
)  


















 # coord_flip()+
  
  # facet_wrap(
  #   vars(!!sym(split_plots_by)),
  #   labeller = labeller(
  #     .default = function(x) paste0(split_plots_clean, ": ", x)
  #   )
  # ) +
  
  # xlim(0,150)+
  #ylim(0,1)+
  
  # scale_color_manual(values = colors_map)+
  
 
  #theme(legend.position = c(0.2, 0.2))
  
  # labs(title=paste0("Overall Survival ", "\n",
  #                   "by: ",  split_curves_clean, "\n",
  #                   "facet: ", split_plots_clean, "\n",
  #                   "subset: ", plot_subset_clean, ": ", plot_subset_values[subset]),
  #      x = "Days Post Injection",
  #      y = "% Survival",
  #      color=split_curves_clean
  #      
  # )+
  # theme(plot.title = element_text(hjust = 0),
  #       aspect.ratio=aspectratio,
  #       plot.margin = margin(5, 5, 40, 5))+




# prop_test_results <- df_props %>%
#   group_by(hist_cat_name, resultant_geno) %>%
#   reframe(
#     p_value = tryCatch({
#       test <- prop.test(x = n, n = total_n)
#       test$p.value
#     }, error = function(e) NA)  # Return NA if there's an error
#   )






