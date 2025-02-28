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


########dataset read in and construction######################
#reading in current dataset. 
coh1<-readRDS("nf1g/surv/cohorts-2025-01-07v2.rds")%>%
  mutate(aod = as.numeric(aod))

df0<-coh1%>%

  #clean up event and aod column. 
  mutate(event = if_else(is.na(exp_end_date),1,0))%>%
  #mutate(aod=if_else(event==0, 150,aod))%>%
  
  #making a convenient variable to easily compare nf1 vs. non-nf1s. 
  mutate(nf1=substr(resultant_geno, 1, 6))%>%
  mutate(geno_bg=substr(resultant_geno, 8, nchar(resultant_geno)))%>%
 
  
  #ensure that if there is a note to exclude from hist, that the hist categorization correctly reflects that.
  mutate(hist=if_else(!is.na(exclude_from_hist_reason)&event==0, "xh-ne", hist))%>%
  mutate(hist=if_else(!is.na(exclude_from_hist_reason)&event==1, "xh-e", hist))%>%
  
  
  #add hist category for not included in hist. categorization
  mutate(hist=if_else(is.na(hist)&event==0, "xh-ne", hist))%>%
  mutate(hist=if_else(is.na(hist)&event==1, "xh-e", hist))%>%
  
  #correcting categorization naming. category 3.1 should be pretumorigenic. 
  #3.15 should be 3.1.
  #correct cat and add metadata indicating change.
  
  #correct 3.1
  mutate(metadata = if_else(hist == "3.1", 
                            if_else(is.na(metadata), "formerly 3.1 graded", 
                                    paste0(metadata, "; formerly 3.1 graded")), 
                            metadata))%>%

  mutate(hist=if_else(hist=="3.1", "PT", hist))%>%
  
  #correct 3.15
  mutate(metadata = if_else(hist == "3.15", 
                            if_else(is.na(metadata), "formerly 3.15 graded", 
                                    paste0(metadata, "; formerly 3.15 graded")), 
                            metadata))%>%
  mutate(hist=if_else(hist=="3.15", "3.1", hist))





df1<-df0%>%
  
  
  #temp selection to vis. df
  #select("mouse_num", "resultant_geno", "hist", "event", "metadata", exclude_from_hist_reason, aod)%>%
  mutate(hist=toupper(hist))%>%
  
  #get 1-letter/number shorthand for histological category.
  mutate(hist_cat=substr(hist, 1,1))%>%

  #get hist_grade from hist information column. note that excluded from histology will
  #still show a "grade" but it's not actually a grade. 
  mutate(hist_grade=substr(hist, 3,4))%>%
  mutate(hist_grade=gsub("\\.", "", hist_grade))%>%
  mutate(hist_catgrade=paste0(hist_cat, ".", hist_grade))%>%
  
  
  mutate(across(c(hist_cat, hist_grade, hist_catgrade), 
               ~ if_else(hist == "NED" | grepl("^X", hist)| grepl("^PT", hist), hist, .x)))%>%
  
  
  mutate(hist_grade=if_else(hist_grade=="N", "NGA", hist_grade))
    


sheet_id<-"https://docs.google.com/spreadsheets/d/1SALlvQdI5zyhZLgUg3XS1kmb2M8Skr_8Pjt9y_spWGE/edit?gid=667454129#gid=667454129"
gs4_auth(email = "gotrimski@gmail.com") 




s_hist<-read_sheet(sheet_id, sheet = "hist") %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~ na_if(.x, "NULL"))) %>%
  janitor::clean_names() %>%
  select(1:3) %>%
  filter(!is.na(.[[3]])) %>%
  mutate(across(3, as.numeric))%>%
  arrange(.[[3]])%>%
  mutate(across(2, ~ factor(.x, levels = .x))) %>% # Preserve order as factor levels
  select(-3)
  
  

s_hist_cat <- read_sheet(sheet_id, sheet = "hist_cat") %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~ na_if(.x, "NULL"))) %>%
  janitor::clean_names() %>%
  select(1:3) %>%
  filter(!is.na(.[[3]])) %>%
  mutate(across(3, as.numeric))%>%
  arrange(.[[3]])%>%
  mutate(across(2, ~ factor(.x, levels = .x))) %>% # Preserve order as factor levels
  select(-3)



s_hist_grade <- read_sheet(sheet_id, sheet = "hist_grade") %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~ na_if(.x, "NULL"))) %>%
  janitor::clean_names() %>%
  select(1:3) %>%
  filter(!is.na(.[[3]])) %>%
  mutate(across(3, as.numeric))%>%
  arrange(.[[3]])%>%
  mutate(across(2, ~ factor(.x, levels = .x))) %>% # Preserve order as factor levels
  select(-3)



s_hist_catgrade<-read_sheet(sheet_id, sheet = "hist_catgrade") %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~ na_if(.x, "NULL"))) %>%
  janitor::clean_names() %>%
  select(1:3) %>%
  filter(!is.na(.[[3]])) %>%
  mutate(across(3, as.numeric))%>%
  arrange(.[[3]])%>%
  mutate(across(2, ~ factor(.x, levels = .x))) %>% # Preserve order as factor levels
  select(-3)






df2<-df1%>%
  left_join(s_hist)%>%
  left_join(s_hist_cat)%>%
  left_join(s_hist_grade)%>%
  left_join(s_hist_catgrade)%>%

  relocate(mouse_num, hist, hist_cat, hist_grade, hist_catgrade,
           hist_name, hist_cat_name, hist_grade_name, hist_catgrade_name,
           exclude, metadata)%>%


  mutate(full_cohort="1")%>%
  
  
  mutate(resultant_geno_name=strain_injection)





saveRDS(df2, "nf1g/surv/cohorts-2025-02-27.rds")












