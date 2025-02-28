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

num_categories<-6
curve_category_names<-readRDS("nf1g/surv/cohorts-2025-02-27.rds")%>%select(resultant_geno)%>%arrange(resultant_geno)%>%unique()%>%pull()


geno_colors_map <- list(resultant_geno=setNames(hue_pal()(num_categories), curve_category_names))



geno_colors_map<-list(resultant_geno=setNames(c("#B79F00",
                                              "#F8766D", 
                                              "#F564E3",
                                              "#619CFF",
                                              "#00BFC4", 
                                              "#00BA38"), curve_category_names))

saveRDS(geno_colors_map, "nf1g/ds/colors_map_geno.rds")








library(googlesheets4)
sheet_id<-"https://docs.google.com/spreadsheets/d/1SALlvQdI5zyhZLgUg3XS1kmb2M8Skr_8Pjt9y_spWGE/edit?gid=667454129#gid=667454129"
gs4_auth(email = "gotrimski@gmail.com") 


for (sn in 1:length(sheet_names(sheet_id))){


name<-sheet_names(sheet_id)[sn]  
  
assign(x=paste0("s_",name),
  
  value=read_sheet(sheet_id, sheet = name) %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~ na_if(.x, "NULL"))) %>%
  janitor::clean_names() %>%
  select(2,5) %>%
  filter(!is.na(.[[1]]))
  
)
}






