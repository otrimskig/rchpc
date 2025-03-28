source("libs.R")
library(dtplyr)
suppressPackageStartupMessages(library(ComplexHeatmap))
ht_opt$message = FALSE
library(gridtext)
library(circlize)
library(purrr)
library(tidyverse)


gsva_u0<-readRDS("nf1g/ds/gsva/gsva_pathways_matrix.rds")

gsva_meta<-readRDS("nf1g/ds/gsva/gsva_pathways_meta.rds")

sample_info<-readRDS("nf1g/ds/v10-per_sample_updated.rds")




patho_cat_names<-sample_info%>%
  select(patho_cat_name)%>%
  unique()%>%
  pull()%>%
  as.character()






p<-3

patho_cat_names[p]

