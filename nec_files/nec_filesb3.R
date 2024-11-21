source("libs.R")

library(tidyverse)
library(readxl)
library(janitor)
library(purrr)
library(lubridate)
library(stringr)


a<-readRDS("nec_files/2024-11-14-combined_cleaned_necs.rds")%>%
 mutate(index=1:n())%>%
relocate(index)%>%


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



b<-readRDS("nec_files/2024-11-14-combined_cleaned_necs-part2.rds")



c1<-bind_rows(a,b)



c2<-c1%>%
   
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





  arrange(desc(necropsy_date))%>%
  mutate(index=1:n())%>%

  relocate(bd, sac_date, .after = necropsy_date)%>%





  unite(
    treatment_duration,
    c("treatment_duration"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  mutate(treatment_duration=if_else(treatment_duration=="", NA, treatment_duration))%>%


  unite(
    virus_1,
    c("virus_1"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  mutate(virus_1=if_else(virus_1=="", NA, virus_1))%>%




  unite(
    duration_of_dox,
    c("duration_of_dox", "duration_of_dox_tx"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  mutate(duration_of_dox=if_else(duration_of_dox=="", NA, duration_of_dox))%>%

  unite(
    blood_collection,
    c("blood_collection"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  mutate(blood_collection=if_else(blood_collection=="", NA, blood_collection))%>%



  unite(
    cell_line_injected,
    c("cell_line_injected"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  mutate(cell_line_injected=if_else(cell_line_injected=="", NA, cell_line_injected))%>%


  

  unite(
    necropsy_date,
    c("necropsy_date", "todays_date"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  mutate(necropsy_date=if_else(necropsy_date=="", NA, necropsy_date))%>%







  
  mutate(dox_tx=if_else(tolower(dox_tx)=="no", NA, dox_tx))%>%


  unite(
    treatment,
    c("treatment"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  mutate(treatment=if_else(treatment=="", NA, treatment))%>%








  
  
  





















  mutate(bd=suppressWarnings(as_date(bd)))%>%
  

  select(where(~ !all(is.na(.))))




saveRDS(c2, "nec_files/2024-11-14-combined_cleaned_necs-part2.rds")



writexl::write_xlsx(c2, "nec_files/necropsy_files-plus-outside.xlsx")







# 
#
# 
# writexl::write_xlsx(x6, "nec_files/necropsy_files2.xlsx")
# 
# 
# 
# 
# 
# #select(sort(names(.)))