library(tidyverse)
library(dtplyr)
library(ComplexHeatmap)
library(gridtext)
library(purrr)



gsva_u<-readRDS("rebecca/gsva_u.rds")

#subset however you like, based on any analysis.
gsva_u_subset<-gsva_u[1:50,]


