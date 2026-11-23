source("libs.R")

library(tidyverse)
library(readxl)
library(janitor)
library(purrr)

x2.1<-readRDS("nec_files/26nec0.rds")%>%
  mutate(across(everything(), as.character))%>%
  mutate(across(everything(), ~na_if(., "NA")))

necdb26<-readRDS("nec_files/2026-06-necdb.rds")%>%
  mutate(across(everything(), as.character))%>%
  mutate(across(everything(), ~na_if(., "NA")))
  




convert_excel_date_column <- function(x) {
  
  x_chr <- trimws(as.character(x))
  
  out <- suppressWarnings(as.Date(
    as.numeric(x_chr),
    origin = "1899-12-30"
  ))
  
  # Fill in failed conversions with parsed text dates
  missing <- is.na(out)
  
  out[missing] <- suppressWarnings(as.Date(
    parse_date_time(
      x_chr[missing],
      orders = c(
        "ymd", "mdy", "dmy",
        "m/d/Y", "d/m/Y",
        "m-d-Y", "d-m-Y"
      )
    )
  ))
  
  out
}


x2.2<-x2.1%>%
  mutate(across(todays_date, convert_excel_date_column))%>%
  mutate(across(bd, convert_excel_date_column))%>%
  mutate(across(sac_date, convert_excel_date_column))%>%
  
  mutate(date_tumor_reported_note=date_tumor_reported)%>%
  mutate(across(date_tumor_reported, convert_excel_date_column))%>%
  mutate(date_tumor_reported_note=if_else(as.character(date_tumor_reported)==date_tumor_reported_note, 
         NA_character_, as.character(date_tumor_reported_note)))%>%
  
  
  select(sort(names(.)))



x2.3<-x2.2%>%
  
  
  unite(
   additional_comments,
    c("additional_comments", "na", starts_with("na_")),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  mutate(additional_comments=if_else(additional_comments=="", NA, additional_comments))%>%
  mutate(across(everything(), as.character))
  
  
 






necdb0<-necdb26%>%
  
  bind_rows(x2.3)%>%
  
  select(sort(names(.)))%>%
  
  
  unite(
    disposition_notes,
    c(starts_with("disposition_notes")),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  
 
  unite(
    age_at_sac,
    c(starts_with("age_at_sac")),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  
  unite(
    by,
    c("by", "by_2"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  
  
  unite(
    date_tumor_reported_note,
    c(starts_with("date_tumor_reported_note")),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  
  
  unite(
    duration_of_dox,
    c(starts_with("duration_of_dox")),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  
  
  unite(
    injection_date,
    c("injection_date", "injection_date_1"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  
  
  unite(
    additional_comments,
    c("additional_comments", "x3"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  
  rename(injected_by=by)%>%
  rename(date_of_today=todays_date)%>%
  rename(date_of_necropsy=necropsy_date)%>%
  rename(pfa_transcardial_perfusion=x4_percent_pfa_transcardial_perfusion)%>%
  
  rename(date_of_birth=bd)%>%
  rename(date_of_death=death_date)%>%
  
  rename(necropsy_by=initials)%>%
  
  
  
  unite(
    date_of_necropsy,
    c("date_of_necropsy", "date_of_today"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  
  unite(
    original_filename,
    c("original_filename", "filename"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  
  
  unite(
    genotype,
    c("genotype", "strain"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  
  
  
  unite(
    date_of_death,
    c("date_of_death", "sac_date"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  
  select(sort(names(.)))%>%
  mutate(across(everything(), as.character))%>%
  mutate(across(everything(), ~na_if(., "")))







saveRDS(necdb0, "nec_files/26necdb.rds")

necdb0.1<-necdb0%>%
  
  rename_with(~ gsub("_", " ", .x))

writexl::write_xlsx(necdb0.1, "nec_files/26necdb.xlsx")


