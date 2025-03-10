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

library(googlesheets4)
gs4_auth(email = "gotrimski@gmail.com")

df0<-read_sheet("https://docs.google.com/spreadsheets/d/1uZOQGLUsoWWLQ0_zj-sHzGbVjIo7NVnROKd4z951Dw8/edit?gid=1780804155#gid=1780804155", 
                sheet = "Sheet3")%>%
  mutate(across(1:last_col(), as.character))%>%
  
  #then replace all "NULL" with NA. 
  mutate(across(1:last_col(), function(x){na_if(x, "NULL")}))%>%
  
  janitor::clean_names()%>%
  
  mutate(aod=as.numeric(aod))%>%
  mutate(event=as.numeric(event))






#use to create factors and determine factor levels. 
source("nf1g/colors_map_create.R")
#check and edit levels (order here: https://docs.google.com/spreadsheets/d/1SALlvQdI5zyhZLgUg3XS1kmb2M8Skr_8Pjt9y_spWGE/edit?gid=763297397#gid=763297397)

col_map<-readRDS("nf1g/surv/colors_map_surv.rds")


df1 <- df0 %>%
  mutate(across(names(col_map), ~ factor(.x, levels = names(col_map[[cur_column()]]))))









saveRDS(df1, "nf1g/surv/cohorts-2025-03-07-v2.rds")
