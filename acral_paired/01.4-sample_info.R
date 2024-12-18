library(tidyverse)
library(janitor)
library(stringr)


a<-read_tsv("acral_paired/ds/gnomex_exp_data.txt")%>%
  select(where(~ !all(is.na(.))))%>%
  clean_names()%>%
  relocate(id, sample_name)%>%
 
  
  
  mutate(sample_num=sub("25622X", "",id))%>%
  
  relocate(sample_num)%>%
  mutate(sample_num=sprintf("%02d", as.numeric(sample_num)))%>%
  
  mutate(sample_id= paste0("x25622x", sample_num))%>%
  
  relocate(sample_id)%>%
  select(-sample_num, -id)%>%
  
  mutate(mouse_num = substr(sample_name, 1, 5))%>%
  mutate(sample_type=substr(sample_name, 7,8))%>%
  mutate(sample_type=if_else(sample_type=="fp", "acral", "subq"))%>%
  
  select(-sample_name)%>%
  relocate(mouse_num, sample_type)
  
