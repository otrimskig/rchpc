source("libs.R")
suppressMessages(source("nf1g/colors_map_create.R", echo = FALSE))

library(tidyverse)
library(ggpubr)
library(gridExtra)
library(ggkm)
library(scales)
library(survival)
library(survminer)
library(RColorBrewer)

source("ggplot_draw_square.R")


########dataset read in and construction######################
#reading in current dataset. 
coh1<-readRDS("nf1g/surv/cohorts-2025-03-07-v2.rds")
col_map<-readRDS("nf1g/surv/colors_map_surv.rds")
aspectratio<-.6

df1<-coh1%>%
  #most conservative exclusion criteria
  filter(is.na(exclude))%>%
  filter(!is.na(include_in_surv))
  #filter(is.na(exclude_hist))%>%
  #filter(as.numeric(hist_grade)>1)%>%
  

ct_geno<-df1%>%
  group_by(resultant_geno)%>%
  
  summarise(count=n())


ct_tt<-df1%>%
  group_by(hist_cat_name)%>%
  
  summarise(count=n())

sum_df1<-df1%>%
  group_by(resultant_geno, hist_cat_name, hist_grade_name)%>%
  
  summarise(n_geno=n())%>%
  
  ungroup()%>%
  
  
  
  arrange(resultant_geno, hist_cat_name, hist_grade_name)
 
  
   
  


comp_1<-sum_df1%>%
  filter(grepl("\\d$", hist_grade_name))%>%
  mutate(across(where(is.factor), droplevels))%>%
  complete(resultant_geno, hist_cat_name, hist_grade_name,
           fill = list(n_geno = 0))


comp_2<-sum_df1%>%
  filter(!grepl("\\d$", hist_grade_name))%>%
  #filter(!grepl("No grade assigned", hist_grade_name))
  mutate(across(where(is.factor), droplevels))%>%

  complete(resultant_geno, hist_cat_name,
           fill = list(n_geno = 0))%>%
  
  mutate(hist_grade_name=if_else(hist_cat_name=="Glioneuronal tumors"|
                                   hist_cat_name=="No histological classification", "No grade assigned", hist_grade_name))%>%
  
  mutate(hist_grade_name=if_else(is.na(hist_grade_name), hist_cat_name, hist_grade_name))






comp_3<-full_join(comp_1, comp_2)%>%
  mutate(hist_grade_name = factor(hist_grade_name, levels = names(col_map$hist_grade_name)))












sum_df2<-comp_3%>%
  ungroup%>%group_by(resultant_geno)%>%
  reframe(resultant_geno,
          hist_cat_name, 
          hist_grade_name,
          n_geno, 
          n_total_geno=sum(n_geno),
          prop_geno=n_geno/n_total_geno,
          se_prop = sqrt((prop_geno * (1 - prop_geno)) / n_total_geno))%>%
  
  group_by(hist_cat_name)%>%
  
          
  reframe(across(everything()), n_total_hist=sum(n_geno),
            prop_total_hist=sum(prop_geno),
          prop_prop_total_hist=prop_geno/prop_total_hist,
          se_prop_scaled=se_prop/prop_total_hist)%>%

  relocate(resultant_geno)%>%
  relocate(n_total_hist, .after=n_total_geno)%>%
  
  
  mutate(hist_cat_name_numeric=as.numeric(hist_cat_name))







cat_levels <- unique(sum_df2$hist_grade_name)
sum_df2$hist_grade_name_numeric <- match(sum_df2$hist_grade_name, cat_levels) * 0.2  # Adjust 0.8 to fine-tune spacing

gg_data<-sum_df2%>%
  filter(grepl("\\d$", hist_grade_name)|
           hist_grade_name=="No grade assigned")%>%
  filter(hist_cat_name!="No histological classification")%>%
  mutate(across(where(is.factor), droplevels))%>%
  
  complete(hist_grade_name, hist_cat_name, resultant_geno, fill=list(n_geno=0))%>%
  mutate(prop_geno_dummy=if_else(n_geno==0, -.5, prop_geno*100))
  

ggplot(gg_data) +
  geom_col(aes(x = hist_grade_name_numeric,
               fill = hist_cat_name,
               y = prop_geno_dummy),
           width = 0.4,
           position = position_dodge(width = .5, preserve="total"),
           key_glyph = draw_square)+
  facet_grid(vars(resultant_geno))
  # 
  # 
  # scale_x_continuous(breaks = unique(gg_data$hist_cat_name_numeric),
  #                    labels = unique(gg_data$hist_cat_name))
  # 
  # 




  
  geom_errorbar(aes(
    x = hist_cat_name_numeric, 
    ymin = ifelse(perc > 0, perc - se * 100, NA), 
    ymax = ifelse(perc > 0, perc + se * 100, NA),
    group = resultant_geno_name
  ),
  position = position_dodge(0.5),
  width = .2,
  linewidth = .01,
  alpha=.5)
  