source("libs.R")

library(tidyverse)
library(readxl)
library(janitor)
library(purrr)
library(lubridate)
library(stringr)


c2<-readRDS("nec_files/2024-11-14-combined_cleaned_necs-part2.rds")

library(googlesheets4)
gs4_auth(email = "gotrimski@gmail.com")

tibble(names(c2))%>%
  range_write("https://docs.google.com/spreadsheets/d/1ZZJa9YusmSC5TeI3btwDuwHlZGg_F3BvLQmUL79Bs2E/edit?gid=0#gid=0",
              sheet = "Sheet1",
              .,
              reformat=FALSE,
              range = "A1")



df1<-read_sheet("https://docs.google.com/spreadsheets/d/1ZZJa9YusmSC5TeI3btwDuwHlZGg_F3BvLQmUL79Bs2E/edit?gid=0#gid=0",
           sheet="Sheet1")%>%
  mutate(across(1:last_col(), as.character))%>%
  
  #then replace all "NULL" with NA. 
  mutate(across(1:last_col(), function(x){na_if(x, "NULL")}))%>%
  
  janitor::clean_names()


df2<-df1%>%
  mutate(clean_index=1:n())%>%
  rename(clean1=cleaned_names)%>%
  
  mutate(clean2=gsub("_", " ", tolower(clean1)))%>%
  arrange(names_c2)%>%
  
  mutate(a_index=1:n())


r_names<-df2%>%
  rename(r_names=names_c2)%>%
  select(r_names)


other_names<-df2%>%
  select(-names_c2)%>%
  mutate(r_names=gsub(" ", "_", clean2))


positions<-full_join(other_names, r_names)
  








