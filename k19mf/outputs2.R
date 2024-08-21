source("libs.R")
library(tidyverse)
library(dtplyr)


d<-readRDS("k19mf/dexps/dexp-tumor_type-subq v. foot padexcluded2s.rds")

write_csv(d, "k19mf/ds/dexp-foot_pads-excluded2ms.csv")


e<-readRDS("k19mf/dexps/dexp-tumor_type-subq v. foot pad.rds")


write_csv(e, "k19mf/ds/dexp-foot_pads.csv")



