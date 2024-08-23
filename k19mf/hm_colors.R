library(tidyverse)



all_sample_info<-readRDS("k19mf/ds/vm-00-sample_info.rds")%>%
  select(genotype, tumor_type)

pivot<-all_sample_info%>%
  mutate(across(1:last_col(), ~as.character(.)))%>%
  pivot_longer(cols=1:last_col(), names_to = "variable", values_to="condition")%>%
  unique()%>%
  filter(variable!="total_reads"&variable!="aod"&variable!="mouse_num"&variable!="sample_id")

library(googlesheets4)

gs4_auth(email = "gotrimski@gmail.com")


stop("check gs output before continuing.")
pivot%>%
  arrange(variable, condition)%>%
range_write("https://docs.google.com/spreadsheets/d/1LIu9mXHMAs1-O95S1kZtGF2sJxBpFJgBWkQRsUOU-kU/edit#gid=0",
                          sheet = "k19mf",
                          .,
                           reformat=TRUE,
                          range = "A1")






