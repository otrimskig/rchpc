source("libs.R")

library(tidyverse)
library(readxl)
library(janitor)
library(purrr)
library(lubridate)
library(stringr)


x4<-readRDS("nec_files/2024-11-14-combined_cleaned_necs.rds")%>%
  mutate(index=1:n())%>%
  relocate(index)
  

x5<-x4%>%
  
  filter(initials!="ATTENTION:")%>%
  
  
  
  
  mutate(necropsy_date=if_else(necropsy_date=="", NA, necropsy_date))%>%
  mutate(necropsy_date=if_else(index=="648", "2023-02-21", necropsy_date))%>%
  mutate(necropsy_date=if_else(index=="649", "2023-02-21", necropsy_date))%>%
  mutate(necropsy_date=if_else(index=="1684", "2014-02-11", necropsy_date))%>%
  mutate(necropsy_date=if_else(index=="3033", "2015-12-05", necropsy_date))%>%
  mutate(necropsy_date=if_else(index=="267", "2021-11-15", necropsy_date))%>%
  mutate(necropsy_date=if_else(index=="2198", "2014-10-28", necropsy_date))%>%
  
  
  mutate(necropsy_date=if_else(necropsy_date=="3.20.14", "2014-03-20", necropsy_date))%>%
  mutate(necropsy_date=if_else(necropsy_date=="4/21/203", "2023-04-21", necropsy_date))%>%
  mutate(necropsy_date=if_else(necropsy_date=="4/20/203", "2023-04-20", necropsy_date))%>%
  
  
  mutate(overridden= if_else(grepl("overri", necropsy_date), "og necropsy sheet overwritten", NA))%>%
  
  mutate(necropsy_date=if_else(index=="5351"|index=="5333", "2021-03-15", necropsy_date))%>%
  
  
  mutate(necropsy_date=as_date(necropsy_date))%>%
  
  arrange(desc(necropsy_date))




x6<-x5%>%
  
  mutate(date_treatment_started=as_date(date_treatment_started))%>%
 
  mutate(date_domed_head_reported=if_else(date_domed_head_reported=="--", NA, date_domed_head_reported))%>%
  mutate(date_domed_head_reported=as_date(date_domed_head_reported))%>%
  

  mutate(injection_date_1 = as_date(str_extract(injection_date, "^\\d{4}-\\d{2}-\\d{2}$")))%>%
  mutate(injection_date=if_else(!is.na(injection_date_1), NA, injection_date))%>%
  rename(injection_date_notes=injection_date)%>%
  mutate(injection_date_2=NA)%>%
  
  mutate(date_tumor_reported_1 = as_date(str_extract(date_tumor_reported, "^\\d{4}-\\d{2}-\\d{2}$")))%>%
  mutate(date_tumor_reported=if_else(!is.na(date_tumor_reported_1), NA, date_tumor_reported))%>%
  rename(date_tumor_reported_notes=date_tumor_reported)%>%
  rename(date_tumor_reported=date_tumor_reported_1)%>%
  
  mutate(sac_date_clean = as_date(str_extract(sac_date, "^\\d{4}-\\d{2}-\\d{2}$")))%>%
  mutate(sac_date=if_else(!is.na(sac_date_clean), NA, sac_date))%>%
  rename(sac_date_notes=sac_date)%>%
  rename(sac_date=sac_date_clean)%>%
  
  
  mutate(bd=as_date(bd))%>%
  
  
  arrange(desc(necropsy_date))%>%
  mutate(index=1:n())%>%
  
  relocate(bd, sac_date, .after = necropsy_date)
  

#select(where(~ !all(is.na(.))))

writexl::write_xlsx(x6, "nec_files/necropsy_files2.xlsx")





#select(sort(names(.)))