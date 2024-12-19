library(tidyverse)
library(janitor)
library(stringr)
library(dplyr)


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
  


feature_stats<-readRDS("acral_paired/ds/featurecount_stats.rds")%>%
  rename(bam_file=sample_id)%>%
  mutate(sample_id=tolower(paste0("x", substr(bam_file, 1,8))))%>%
  relocate(sample_id)%>%
  dplyr::select(where(~ !all(. == 0)))
  



b<-left_join(a, feature_stats)

saveRDS(b, "acral_paired/ds/v00-sample_info.rds")
