source("libs.R")
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(grid)
library(ggkm)
library(scales)
library(survival)
library(survminer)
library(RColorBrewer)





num_categories<-6
curve_category_names<-readRDS("nf1g/surv/cohorts-2025-01-07.rds")%>%select(resultant_geno)%>%arrange(resultant_geno)%>%unique()%>%pull()


geno_colors_map <- list(resultant_geno=setNames(hue_pal()(num_categories), curve_category_names))



geno_colors_map<-list(resultant_geno=setNames(c("#B79F00",
                                              "#F8766D", 
                                              "#F564E3",
                                              "#619CFF",
                                              "#00BFC4", 
                                              "#00BA38"), curve_category_names))

saveRDS(geno_colors_map, "nf1g/ds/colors_map_geno.rds")





