library(tidyverse)
library(googlesheets4)

#edit colors here: https://docs.google.com/spreadsheets/d/1LIu9mXHMAs1-O95S1kZtGF2sJxBpFJgBWkQRsUOU-kU/edit#gid=0



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


#read in updated metadata for labels.(this will ensure that this is the 
#most updated version of the metadata). If this needs to be updated,
#use the updating_category_names.R which takes input from ds/nf1-renaming dataset.xlsx and
#outputs to the below file.

updated_meta<-readRDS("ds/v10-per_sample_updated.rds")%>%
  select(starts_with("patho"), starts_with("resultant"))


ud_a<-updated_meta[1:2]%>%unique()
ud_b<-updated_meta[3:4]%>%unique()
ud_c<-updated_meta[5:6]%>%unique()
ud_d<-updated_meta[7:8]%>%unique()


colors<-sheets_df%>%
  pivot_wider(names_from = variable, values_from=condition)%>%
  left_join(ud_a)%>%
  left_join(ud_b)%>%
  left_join(ud_c)%>%
  left_join(ud_d)%>%
  pivot_longer(cols=2:last_col(), names_to="variable", values_to="condition")%>%
  filter(!is.na(condition))


variable_names<-colors%>%
  filter(!grepl("^aod", variable))%>%
  group_by(variable)%>%slice(1)%>%ungroup()%>%select(variable)%>%pull()



colors_list<-list()
for (i in 1:length(variable_names)){
  
  
  v<-colors%>%filter(variable==variable_names[i])
  
  colors_list[[variable_names[i]]]<-setNames(v$hex, v$condition)
  
}



saveRDS(colors_list, "ds/colors_list.rds")






