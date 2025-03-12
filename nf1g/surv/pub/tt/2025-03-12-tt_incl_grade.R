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
  filter(!is.na(include_in_surv))%>%
  filter(is.na(exclude_hist))%>%
  filter(as.numeric(hist_grade)>1)%>%
  select(mouse_num, resultant_geno, hist_cat_name, hist_grade_name)

df1_props<-df1%>%
  group_by(resultant_geno, hist_cat_name)
