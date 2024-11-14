source("libs.R")

library(tidyverse)
library(readxl)
library(janitor)
library(purrr)

library(foreach)

if (!exists("n.cores")) {
  
  "initilizing cores..."
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  doParallel::registerDoParallel(cl = my.cluster)
  
  #check if it is registered (optional)
  foreach::getDoParRegistered()
  
  "parallel cores initialized."
  
}




files<-fs::dir_ls("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/nec_files/")%>%
  as_tibble()


files2<-files%>%
  rename(path=value)%>%
  mutate(fn=basename(path))%>%
  mutate(ft=if_else(grepl("xlsx$", fn), "xlsx", 
                    if_else(grepl("xls$", fn), "xls", NA)))



all_files_list<-list()



x<-foreach(i=1:length(files2$path)) %dopar% {
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(tidyr)

# for(i in 1:length(files2$path)){

read_excel(files2$path[i])%>%
  clean_names()%>%
  mutate(across(everything(), ~ ifelse(. == "`", NA, .)))%>%
  select(where(~ !all(is.na(.))))%>%
  t()%>%
  as_tibble()%>%
  row_to_names(1)%>%
  clean_names()%>%
  
  mutate(filename=files2$fn[i])


}


x1 <- keep(x, ~ nrow(.) >= 1)


x2 <- bind_rows(x1)


x3<-x2%>%select(where(~ !all(is.na(.))))%>%
  
  relocate(filename)%>%
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
  
  
  rename(bx=x)%>%
  
  mutate(bd=if_else(is.na(bd), bx, bd))%>%

  select(-bx)%>%
  
  
  
  unite(
    treatment_duration,
    c("treatment_duration", "treatment_length"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  mutate(treatment_duration=if_else(treatment_duration=="", NA, treatment_duration))%>%


  unite(
    virus_1,
    c("virus_1", "virus1"),
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
    c("blood_collection", "blood_draw"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  mutate(blood_collection=if_else(blood_collection=="", NA, blood_collection))%>%
  
  
  
  unite(
    cell_line_injected,
    c("cell_line_injected", "cells"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  mutate(cell_line_injected=if_else(cell_line_injected=="", NA, cell_line_injected))%>%
  
  
  unite(
   date_of_first_treatment,
    c("date_of_first_treatment", "date_of_first_treatrment"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  mutate(date_of_first_treatment=if_else(date_of_first_treatment=="", NA, date_of_first_treatment))%>%
  
  
  unite(
    necropsy_date,
    c("necropsy_date", "todays_date"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  mutate(date_of_first_treatment=if_else(date_of_first_treatment=="", NA, date_of_first_treatment))%>%
  
  
  mutate(across(everything(), ~ case_when(
    str_to_lower(.) %in% c("na", "n/a") ~ NA_character_,  # Replace with NA if the value matches
    TRUE ~ .  # Keep other values unchanged
  )))%>%
  
  
  select(-ihc_stain_2)%>%
  
  
  
  mutate(treated=if_else(tolower(treated)=="no", NA, treated))%>%
  mutate(dox_tx=if_else(tolower(dox_tx)=="no", NA, dox_tx))%>%
  
  
  unite(
    treatment,
    c("treatment", "treated", "drug_treatment"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  mutate(treatment=if_else(treatment=="", NA, treatment))%>%
  
  unite(
    date_treatment_started,
    c("date_of_first_treatment", "start_of_treatment"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  mutate(date_treatment_started=if_else(date_treatment_started=="", NA, date_treatment_started))%>%
  
  

  
select(where(~ !all(is.na(.))))




saveRDS(x4, "nec_files/2024-11-14-combined_cleaned_necs.rds")



# writexl::write_xlsx(x4, "nec_files/necropsy_files.xlsx")
