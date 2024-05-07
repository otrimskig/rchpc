library(tidyverse)
library(edgeR)



df<-readRDS("23908R/v03-rpkms.rds")%>%
  rename(sample_id=sample_name)

sample_info<-readRDS("23908R/sample_mouse_num.rds")%>%
  mutate(tumor_grouping1 = if_else(is.na(tumor_grouping1), "NED", tumor_grouping1))%>%
  mutate(aod = if_else(is.na(aod), 150, aod))








ms_hist<-read_csv("mouse slide samples - with histology notes - updated_slides.csv")%>%
  mutate(mouse_num = as.character(mouse_num))%>%
  semi_join(sample_info, by="mouse_num")%>%
  
  janitor::clean_names()%>%
  
  filter(!is.na(category_1))

ms_hist2<-ms_hist%>%
  filter(grepl("^\\d", category_1))%>%
  
  select(mouse_num, category_1)%>%
  
  rename(patho_cat_det=category_1)%>%
  
  mutate(patho_cat = substr(patho_cat_det, 1,3))




ms_info<-left_join(sample_info, ms_hist2)%>%
  mutate(patho_cat_det = if_else(is.na(patho_cat_det), "NED", patho_cat_det))%>%
  mutate(patho_cat = if_else(is.na(patho_cat), "NED", patho_cat))



all<-df%>%
  left_join(ms_info, by="sample_id")
  
saveRDS(all, "23908R/v05-all_counts_plus_info.rds")




