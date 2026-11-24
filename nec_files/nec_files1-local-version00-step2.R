source("../libs.R")

library(tidyverse)
library(readxl)
library(janitor)
library(purrr)



x4<-readRDS("2026-11-23-combined_cleaned_necs.rds")

la<-read_excel("26necdb_1_from-LA.xlsx")

la2<-la%>%
  janitor::clean_names()%>%
  select(mouse_id)%>%
  unique()
