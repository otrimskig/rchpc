source("libs.R")
library(tidyverse)
library(dtplyr)
library(ggplot2)


#get vector of all .rds files containing dexp analyses.
all_dexps<-list.files("nf1g/dexps/atrx_comps", full.names = T, pattern="^dexp.*\\.rds$")


readRDS("")


df0<-readRDS(all_dexps[1])




