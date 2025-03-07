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


########dataset read in and construction######################
#reading in current dataset. 
coh1<-readRDS("nf1g/surv/cohorts-2025-02-27.rds")%>%
  mutate(aod = as.numeric(aod))%>%
  mutate(aod=if_else(event==0, 150,aod))



df1<-coh1%>%
  filter(is.na(exclude)|exclude>3)%>%

 filter(is.na(exclude_from_hist_reason)|
          grepl("^fd$", exclude_from_hist_reason)|
          grepl("^fd ", exclude_from_hist_reason)|
          grepl("delay", exclude_from_hist_reason)|
          grepl("^fd$", g_macro_obs)|
          grepl("not collected", exclude_from_hist_reason))%>%
  
  mutate(exclude_hist=if_else(!is.na(exclude_from_hist_reason), "1", NA_character_))%>%
  mutate(include_in_surv="1")%>%
  
  relocate(exclude_hist, exclude, metadata)














te<-coh1%>%
  filter(is.na(exclude)|exclude>3)%>%
  relocate(aod, .after="resultant_geno")%>%
  anti_join(df1, by="mouse_num")%>%
  arrange(resultant_geno, aod)%>%
  
  mutate(exclude=if_else(mouse_num=="25499"|mouse_num=="25500"|mouse_num=="25502", "2", exclude))%>%
  mutate(metadata=if_else(mouse_num=="25499"|mouse_num=="25500"|mouse_num=="25502", 
                         
                           if_else(is.na(metadata), #handle positive results of above in 2 different ways
                                   "non tumor-related reason for death", paste0(metadata, "; non tumor-related reason for death")),
                          
                          metadata)
         )%>%
  
  
  
  mutate(metadata=if_else(grepl("cell line", exclude_from_hist_reason), 
                        
                        if_else(is.na(metadata), #handle positive results of above in 2 different ways
                                "samples used for cell lines", paste0(metadata, "; samples used for cell lines")),
                        
                        metadata))%>%
  
  mutate(exclude=if_else((is.na(exclude)|as.numeric(exclude)>3)&grepl("used for cell lines", metadata), "7", exclude))




  
  # 
  # mutate(metadata=if_else(is.na(metadata), "poor data handling. conflicting data sources or data missing.", metadata))%>%
  # 
  # mutate(exclude=if_else(is.na(exclude)|as.numeric(exclude>3), "3", exclude))



write_csv(te, "nf1g/surv/man_exclusion_check0.csv")



te2<-read_csv("nf1g/surv/man_exclusion_check1 - man_exclusion_check1.csv")%>%
  janitor::clean_names()%>%
  select(mouse_num, exclude, exclude_hist, metadata)%>%
  mutate_at(vars(1:last_col()), ~as.character(.))





coh2<-coh1%>%
  mutate(include_in_surv=NA_character_)%>%
  mutate(exclude_hist=NA_character_)%>%
  rows_update(df1, by="mouse_num")%>%
  rows_update(te2, by="mouse_num")%>%
  relocate(mouse_num, include_in_surv, exclude, exclude_hist, metadata, exclude_from_hist_reason)%>%
  mutate(exclude=if_else(exclude=="5", NA_character_, exclude))%>%
  mutate(include_in_surv=if_else(is.na(exclude)|exclude_hist=="1", "1", include_in_surv))
  


saveRDS(coh2, "nf1g/surv/cohorts-2025-03-06.rds")





