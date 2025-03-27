source("libs.R")
library(tidyverse)
library(dtplyr)
library(purrr)
standard_error <- function(x) {
  sd(x) / sqrt(length(x))
}

#get in metadata for samples.
all_sample_info<-readRDS("nf1g/ds/v10-per_sample_updated.rds")%>%
  filter(coh==1|coh==2)

#filter for cohorts we're after.
tumor_types<-all_sample_info%>%
  filter(patho_cat!="4")%>%
  select(patho_cat_name)%>%
  unique()%>%
  pull()%>%
  as.character()

sample_info<-all_sample_info%>%
  filter(patho_cat_name==tum)





