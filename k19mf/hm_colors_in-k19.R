library(tidyverse)
library(googlesheets4)

#authorize user.
gs4_auth(email = "gotrimski@gmail.com")

#read input sheet
colors_df<-read_sheet("https://docs.google.com/spreadsheets/d/1LIu9mXHMAs1-O95S1kZtGF2sJxBpFJgBWkQRsUOU-kU/edit#gid=0", 
                      sheet = "k19mf")%>%
  mutate(across(1:last_col(), as.character))%>%
  
  #then replace all "NULL" with NA. 
  mutate(across(1:last_col(), function(x){na_if(x, "NULL")}))%>%
  
  janitor::clean_names()%>%
  
  filter(!is.na(variable))%>%
  select(variable, condition, hex)

saveRDS(colors_df, "k19mf/ds/colors_df.rds")


variable_names<-unique(colors_df$variable)


colors_list<-list()
for (i in 1:length(variable_names)){
  
  
  v<-colors_df%>%filter(variable==variable_names[i])
  
  colors_list[[variable_names[i]]]<-setNames(v$hex, v$condition)
  
}


saveRDS(colors_list, "k19mf/ds/colors_list.rds")

