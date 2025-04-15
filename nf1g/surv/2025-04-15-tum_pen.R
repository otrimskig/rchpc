source("libs.R")
library(ggplot2)

coh0<-readRDS("nf1g/surv/cohorts-2025-03-07-v2.rds")


coh1<-coh0%>%
  filter(is.na(exclude))%>%
  group_by(resultant_geno)
