source("libs.R")

library(tidyverse)
library(readxl)
library(janitor)
library(purrr)
library(lubridate)
library(stringr)

x<-readRDS("nec_files/outside_files.rds")





x1r <- keep(x, ~ nrow(.) == 1)


x2r <- keep(x, ~ nrow(.) == 2)

x3r<-keep(x, ~ nrow(.) ==3 & length(.)<50)

x2r<-c(x2r, x3r)




xcrap<-keep(x, ~ nrow(.) > 3)






x1r.c <- map(x1r, ~ mutate(.x, across(everything(), as.character)))

for (i in 1:length(x1r.c)){

x1r.c[[i]]<-x1r.c[[i]]%>%
  
  mutate(initials_added=names(x1r.c[[i]])[1])%>%
  select(where(~ !all(is.na(.))))


}


x1.b <- bind_rows(x1r.c)









for (i in 1:length(x2r)){

x2r[[i]]<-x2r[[1]]%>%
  
  mutate(across(everything(), ~ ifelse(. == "`", NA, .)))%>%
  slice(1)%>%
  select(where(~ !all(is.na(.))))



}


x2r.c <- map(x2r, ~ mutate(.x, across(everything(), as.character)))

x2.b <- bind_rows(x2r.c)







x2<-bind_rows(x1.b,  x2.b)










x3<-x2%>%
  
  janitor::clean_names()%>%
  
  mutate(across(everything(), ~ case_when(
    str_to_lower(.) %in% c("na", "n/a") ~ NA_character_,  # Replace with NA if the value matches
    TRUE ~ .  # Keep other values unchanged
  )))%>%
  
  
  
  
  select(where(~ !all(is.na(.))))%>%
  
  


  rename(sheet_name=sheet_name2)%>%
  relocate(filename, sheet_name)%>%
  unite(
    disposition_notes,
    starts_with("dispos"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  mutate(disposition_notes=if_else(disposition_notes=="", NA, disposition_notes))%>%
  
  unite(
    additional_comments,
    starts_with("additional_com"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  mutate(additional_comments=if_else(additional_comments=="", NA, additional_comments))%>%
  
  
  mutate(across(everything(), ~ {
    # Apply logic to each value in the column
    sapply(., function(x) {
      # Check if the value is a 5-digit number between 40000 and 46000
      if (grepl("^\\d{5}$", as.character(x)) && as.numeric(x) >= 40000 && as.numeric(x) <= 46000) {
        # Convert the value into an Excel date and then back to character
        return(as.character(as.Date(as.numeric(x), origin = "1899-12-30")))
      }
      # Return the value as-is if it doesn't meet the condition
      return(as.character(x))
    })
  }))%>%
  
  
  
  mutate(genotype=if_else(is.na(genotype)&grepl("^Nes", na),  na, genotype))%>%
  mutate(na=if_else(genotype==na, NA, genotype))






x4<-x3%>%
  
  
  unite(
    additional_comments,
    c("additional_comments", "na", starts_with("na_")),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  mutate(additional_comments=if_else(additional_comments=="", NA, additional_comments))%>%
  
  
  
  
  mutate(initials=if_else(initials=="initials", NA, initials))%>%
  
  
  
  unite(
    initials,
    c("initials", "initials_added"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  mutate(initials=if_else(initials=="", NA, initials))%>%
  
  
  
  
  
  # # rename(bx=x)%>%
  # 
  # mutate(bd=if_else(is.na(bd), bx, bd))%>%
  # 
  # select(-bx)%>%
  
  
  
  # unite(
  #   treatment_duration,
  #   c("treatment_duration", "treatment_length"),
  #   sep = "; ",
  #   na.rm = TRUE,
  #   remove = TRUE
  # )%>%
  # mutate(treatment_duration=if_else(treatment_duration=="", NA, treatment_duration))%>%
  
  
  # unite(
  #   virus_1,
  #   c("virus_1", "virus1"),
  #   sep = "; ",
  #   na.rm = TRUE,
  #   remove = TRUE
  # )%>%
  # mutate(virus_1=if_else(virus_1=="", NA, virus_1))%>%
  
  
  
  
  # unite(
  #   duration_of_dox,
  #   c("duration_of_dox", "duration_of_dox_tx"),
  #   sep = "; ",
  #   na.rm = TRUE,
  #   remove = TRUE
  # )%>%
  # mutate(duration_of_dox=if_else(duration_of_dox=="", NA, duration_of_dox))%>%
  
  # unite(
  #   blood_collection,
  #   c("blood_collection", "blood_draw"),
  #   sep = "; ",
  #   na.rm = TRUE,
  #   remove = TRUE
  # )%>%
  # mutate(blood_collection=if_else(blood_collection=="", NA, blood_collection))%>%
  # 
  
  
  # unite(
  #   cell_line_injected,
  #   c("cell_line_injected", "cells"),
  #   sep = "; ",
  #   na.rm = TRUE,
  #   remove = TRUE
  # )%>%
  # mutate(cell_line_injected=if_else(cell_line_injected=="", NA, cell_line_injected))%>%
  # 
  
  # unite(
  #   date_of_first_treatment,
  #   c("date_of_first_treatment", "date_of_first_treatrment"),
  #   sep = "; ",
  #   na.rm = TRUE,
  #   remove = TRUE
  # )%>%
  # mutate(date_of_first_treatment=if_else(date_of_first_treatment=="", NA, date_of_first_treatment))%>%
  # 
  
  # unite(
  #   necropsy_date,
  #   c("necropsy_date", "todays_date"),
  #   sep = "; ",
  #   na.rm = TRUE,
  #   remove = TRUE
  # )%>%
  # mutate(date_of_first_treatment=if_else(date_of_first_treatment=="", NA, date_of_first_treatment))%>%
  
  
  

  
  
  
  # mutate(treated=if_else(tolower(treated)=="no", NA, treated))%>%
  # mutate(dox_tx=if_else(tolower(dox_tx)=="no", NA, dox_tx))%>%
  # 
  # 
  # unite(
  #   treatment,
  #   c("treatment", "treated", "drug_treatment"),
  #   sep = "; ",
  #   na.rm = TRUE,
  #   remove = TRUE
  # )%>%
  # mutate(treatment=if_else(treatment=="", NA, treatment))%>%
  # 
  # unite(
  #   date_treatment_started,
  #   c("date_of_first_treatment", "start_of_treatment"),
  #   sep = "; ",
  #   na.rm = TRUE,
  #   remove = TRUE
  # )%>%
  # mutate(date_treatment_started=if_else(date_treatment_started=="", NA, date_treatment_started))%>%
  # 
  # 
  
  
  select(where(~ !all(is.na(.))))




saveRDS(x4, "nec_files/2024-11-14-combined_cleaned_necs-part2.rds")



# writexl::write_xlsx(x4, "nec_files/necropsy_files.xlsx")













# 
# x4<-readRDS("nec_files/2024-11-14-combined_cleaned_necs.rds")%>%
#   mutate(index=1:n())%>%
#   relocate(index)
#   
# 
# x5<-x4%>%
#   
#   filter(initials!="ATTENTION:")%>%
#   
#   
#   
#   
#   mutate(necropsy_date=if_else(necropsy_date=="", NA, necropsy_date))%>%
#   mutate(necropsy_date=if_else(index=="648", "2023-02-21", necropsy_date))%>%
#   mutate(necropsy_date=if_else(index=="649", "2023-02-21", necropsy_date))%>%
#   mutate(necropsy_date=if_else(index=="1684", "2014-02-11", necropsy_date))%>%
#   mutate(necropsy_date=if_else(index=="3033", "2015-12-05", necropsy_date))%>%
#   mutate(necropsy_date=if_else(index=="267", "2021-11-15", necropsy_date))%>%
#   mutate(necropsy_date=if_else(index=="2198", "2014-10-28", necropsy_date))%>%
#   
#   
#   mutate(necropsy_date=if_else(necropsy_date=="3.20.14", "2014-03-20", necropsy_date))%>%
#   mutate(necropsy_date=if_else(necropsy_date=="4/21/203", "2023-04-21", necropsy_date))%>%
#   mutate(necropsy_date=if_else(necropsy_date=="4/20/203", "2023-04-20", necropsy_date))%>%
#   
#   
#   mutate(overridden= if_else(grepl("overri", necropsy_date), "og necropsy sheet overwritten", NA))%>%
#   
#   mutate(necropsy_date=if_else(index=="5351"|index=="5333", "2021-03-15", necropsy_date))%>%
#   
#   
#   mutate(necropsy_date=as_date(necropsy_date))%>%
#   
#   arrange(desc(necropsy_date))
# 
# 
# 
# 
# x6<-x5%>%
#   
#   mutate(date_treatment_started=as_date(date_treatment_started))%>%
#  
#   mutate(date_domed_head_reported=if_else(date_domed_head_reported=="--", NA, date_domed_head_reported))%>%
#   mutate(date_domed_head_reported=as_date(date_domed_head_reported))%>%
#   
# 
#   mutate(injection_date_1 = as_date(str_extract(injection_date, "^\\d{4}-\\d{2}-\\d{2}$")))%>%
#   mutate(injection_date=if_else(!is.na(injection_date_1), NA, injection_date))%>%
#   rename(injection_date_notes=injection_date)%>%
#   mutate(injection_date_2=NA)%>%
#   
#   mutate(date_tumor_reported_1 = as_date(str_extract(date_tumor_reported, "^\\d{4}-\\d{2}-\\d{2}$")))%>%
#   mutate(date_tumor_reported=if_else(!is.na(date_tumor_reported_1), NA, date_tumor_reported))%>%
#   rename(date_tumor_reported_notes=date_tumor_reported)%>%
#   rename(date_tumor_reported=date_tumor_reported_1)%>%
#   
#   mutate(sac_date_clean = as_date(str_extract(sac_date, "^\\d{4}-\\d{2}-\\d{2}$")))%>%
#   mutate(sac_date=if_else(!is.na(sac_date_clean), NA, sac_date))%>%
#   rename(sac_date_notes=sac_date)%>%
#   rename(sac_date=sac_date_clean)%>%
#   
#   
#   mutate(bd=as_date(bd))%>%
#   
#   
#   arrange(desc(necropsy_date))%>%
#   mutate(index=1:n())%>%
#   
#   relocate(bd, sac_date, .after = necropsy_date)
#   
# 
# #select(where(~ !all(is.na(.))))
# 
# writexl::write_xlsx(x6, "nec_files/necropsy_files2.xlsx")
# 
# 
# 
# 
# 
# #select(sort(names(.)))