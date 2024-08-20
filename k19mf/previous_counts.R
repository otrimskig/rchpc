source("libs.R")

library(tidyverse)
library(dtplyr)

setwd(here::here())

k1<-data.table::fread("k19mf/ds/GSE122781_star_counts22.txt")

k2<-data.table::fread("k19mf/ds/GSE122781_series_matrix.txt", fill=TRUE)



k3<-k1%>%
  select(-1)

k4<-summarise_all(k3, sum, na.rm = TRUE)%>%t()%>%
  data.frame()%>%
  janitor::clean_names()%>%
  rename(total_counts_sample=x)%>%
  rownames_to_column("sample")%>%
  mutate(sample=paste0("x", tolower(sample)))%>%
  mutate(sample=if_else(nchar(sample)==8, paste0(substr(sample, 1,7), "0", substr(sample, 8, 8)), sample))


saveRDS(k4, "k19mf/ds/dk_alignment_stats_bowtie.rds")

         