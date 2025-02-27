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


#most conservative exclusion criteria
filter(is.na(exclude)|exclude>3)




















df0<-coh1%>%

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



nf1_renaming_dataset1 <- readxl::read_excel("nf1g/surv/nf1-renaming dataset.xlsx")
nf1_renaming_dataset2 <- readxl::read_excel("nf1g/surv/nf1-renaming dataset.xlsx", sheet = "resultant_geno")
nf1_renaming_dataset3 <- readxl::read_excel("nf1g/surv/nf1-renaming dataset.xlsx", sheet = "patho_cat")
#nf1_renaming_dataset4 <- readxl::read_excel("nf1g/surv/nf1-renaming dataset.xlsx", sheet = "patho_cat2")
nf1_renaming_dataset5 <- readxl::read_excel("nf1g/surv/nf1-renaming dataset.xlsx", sheet = "patho_cat_det")%>%
  rename(hist=patho_cat_det)

df1<-df0%>%
  left_join(nf1_renaming_dataset1)%>%
  left_join(nf1_renaming_dataset2)%>%
  left_join(nf1_renaming_dataset3)%>%
  #left_join(nf1_renaming_dataset4)
  left_join(nf1_renaming_dataset5)%>%
  mutate(full_cohort="1")%>%
  mutate(patho_cat_name=if_else(hist=="3.1", "Pretumorigenic", patho_cat_name))%>%



mutate(hist=if_else(!is.na(exclude_from_hist_reason)&event==1, "excluded from hist (event)", hist))%>%
  mutate(hist=if_else(!is.na(exclude_from_hist_reason)&event==0, "excluded from hist (no event)", hist))%>%



name<-df1%>%select(patho_cat_name)%>%unique()


name_order<-c("Gliomas", "Nerve sheath tumors","Spindle and epithelioid tumors",
"Glioneuronal tumors",    "Pretumorigenic",  "No histological classification", "No evidence of disease")       




df1<-df1%>%
  mutate(patho_cat_name=factor(patho_cat_name, levels=name_order))



aspectratio<-.6















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


