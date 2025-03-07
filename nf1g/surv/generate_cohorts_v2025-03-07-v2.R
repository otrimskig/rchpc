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




saveRDS(df0, "nf1g/surv/cohorts-2025-03-07-v2.rds")
