source("libs.R")

library(tidyverse)
library(readxl)
library(fs)
library(dtplyr)
# library(janitor)
# library(purrr)
# library(lubridate)
# library(stringr)

get_file_extension <- function(filepath) {
  if (grepl("\\.", filepath)) {
    return(sub(".*\\.", "", filepath))
  } else {
    return("")
  }
}
xd<-read_csv("nec_files/x_paths.csv")

xd2<-xd%>%
  janitor::clean_names()%>%
  rename(len=length)%>%
  
  mutate(path=gsub("\\\\", "/", path))

xd3<-xd2%>%  

  mutate(fn=basename(path))%>%
  
  mutate(ft=sapply(fn, get_file_extension))
  
  


saveRDS(xd3, "nec_files/x_index20241114.rds")





se<-xd3%>%


  filter(grepl("necrop", tolower(fn)))%>%
  
  filter(grepl(".xlsx$", fn)|grepl(".xls$", fn))%>%
  
  
  
  filter(!grepl("^X:/Holmen Lab/Necropsy files/", path))



saveRDS(se, "nec_files/files_list.rds")
