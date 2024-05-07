setwd("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig")

library(tidyverse)
library(googlesheets4)

#authorize user.
gs4_auth(email = "gotrimski@gmail.com")

#read input sheet
sheets_df<-read_sheet("https://docs.google.com/spreadsheets/d/1LIu9mXHMAs1-O95S1kZtGF2sJxBpFJgBWkQRsUOU-kU/edit#gid=0", 
                      sheet = "Sheet1")%>%
  mutate(across(1:last_col(), as.character))%>%
  
  #then replace all "NULL" with NA. 
  mutate(across(1:last_col(), function(x){na_if(x, "NULL")}))%>%
  
  janitor::clean_names()%>%
  
  filter(!is.na(variable))%>%
  select(variable, condition, hex)

saveRDS(sheets_df, "colors.rds")



