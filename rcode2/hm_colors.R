setwd("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig")

library(tidyverse)



all_sample_info<-readRDS("23908R/v07-per_sample_info.rds")

pivot<-all_sample_info%>%
  mutate(across(1:last_col(), ~as.character(.)))%>%
  pivot_longer(cols=1:last_col(), names_to = "variable", values_to="condition")%>%
  unique()%>%
  filter(variable!="total_reads"&variable!="aod"&variable!="mouse_num"&variable!="sample_id")

library(googlesheets4)

gs4_auth(email = "gotrimski@gmail.com")

pivot%>%
  arrange(variable, condition)%>%
range_write("https://docs.google.com/spreadsheets/d/1LIu9mXHMAs1-O95S1kZtGF2sJxBpFJgBWkQRsUOU-kU/edit#gid=0",
                          sheet = "Sheet1",
                          .,
                           reformat=TRUE,
                          range = "A1")






