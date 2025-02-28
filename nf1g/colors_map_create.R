source("libs.R")
library(tidyverse)
library(RColorBrewer)


library(googlesheets4)
sheet_id<-"https://docs.google.com/spreadsheets/d/1SALlvQdI5zyhZLgUg3XS1kmb2M8Skr_8Pjt9y_spWGE/edit?gid=667454129#gid=667454129"
gs4_auth(email = "gotrimski@gmail.com") 

col_list<-list()
for (sn in 1:length(sheet_names(sheet_id))){


name<-sheet_names(sheet_id)[sn]  
  
col_df<-read_sheet(sheet_id, sheet = name) %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~ na_if(.x, "NULL"))) %>%
  janitor::clean_names() %>%
  select(2,5) %>%
  filter(!is.na(.[[1]]))%>%  
  rename(hex_col = colnames(.)[2]) 
  
col_list[[colnames(col_df)[1]]] <- setNames(col_df$hex_col, col_df[[colnames(col_df)[1]]])


}



#also save a list of the same colors as resulant_geno, for convenience of conversion.
name<-"resultant_geno"  
col_df<-read_sheet(sheet_id, sheet = name) %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~ na_if(.x, "NULL"))) %>%
  janitor::clean_names() %>%
  select(1,5) %>%
  filter(!is.na(.[[1]]))%>%  
  rename(hex_col = colnames(.)[2]) 
col_list[[colnames(col_df)[1]]] <- setNames(col_df$hex_col, col_df[[colnames(col_df)[1]]])




saveRDS(col_list, "nf1g/surv/colors_map_surv.rds")
