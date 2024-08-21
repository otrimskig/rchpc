source("libs.R")
library(tidyverse)
library(dtplyr)
library(ggplot2)


dexp<-readRDS("k19mf/dexps/dexp-tumor_type-subq v. foot pad.rds")

dexpf<-dexp%>%
  filter(FDR<0.000001)%>%
  filter(logCPM>3)


ggplot(dexpf,aes(x=logFC, y=gene_name_ms))+
  geom_()
