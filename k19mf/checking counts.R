source("libs.R")

library(tidyverse)
library(dtplyr)

setwd(here::here())


m<-readRDS("k19mf/ds/vm-00-all_counts_sample_names_fixed.rds")%>%
  select(-1, -2)%>%
  summarise_all(sum, na.rm=TRUE)%>%t()%>%
  data.frame()%>%janitor::clean_names()%>%
  
  rownames_to_column()%>%
  rename(sample=rowname, total_counts_sample_star=x)


old<-readRDS("k19mf/ds/dk_alignment_stats_bowtie.rds")

comp<-left_join(m, old)


