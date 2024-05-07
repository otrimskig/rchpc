library(fs)
library(tidyverse)

tb<-dir_info(path="normalized_cancer_sets/", recurse = TRUE, type="any")%>%
  as_tibble()

fi<-tb%>%
  mutate(bn=basename(path))%>%
  filter(grepl("txt$", bn))%>%
  
  
  rename(path1=path)%>%
  
  mutate(path2= paste0("normalized_cancer_sets/", bn))%>%
  select(path1,path2)

file_move(fi$path1, fi$path2)

