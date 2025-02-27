#R input using Data to R sheet
#necessary libraries
library(tidyverse)
library(googlesheets4)


##################################
#sheet id and sheet name to read in from sheets. 

sheet_id<-"https://docs.google.com/spreadsheets/d/1RYtJR9CDKXw0JWHGN1nWcqCCIOBYPpcNaew7DsWGs6U/edit?gid=278364251#gid=278364251"
name_of_sheet<-"cohorts_clean_hist"

#authorize user.
gs4_auth(email = "gotrimski@gmail.com")

#read input sheet
sheets_df<-read_sheet(sheet_id, 
                      sheet = name_of_sheet)%>%
  mutate(across(1:last_col(), as.character))%>%
  
  #then replace all "NULL" with NA. 
  mutate(across(1:last_col(), function(x){na_if(x, "NULL")}))%>%
  
  janitor::clean_names()




sheets_df%>%
  group_by(mouse_num)%>%count()%>%
  filter(n>1)



coh1<-sheets_df

saveRDS(coh1, "nf1g/surv/cohorts-2025-01-07v2.rds")








