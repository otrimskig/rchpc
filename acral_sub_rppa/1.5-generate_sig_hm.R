source("libs.R")
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)




mat0<-readRDS("acral_sub_rppa/ds/umat-per_sample_375_rel_diff_th.rds")
stats0<-readRDS("acral_sub_rppa/ds/stats_fr_threshold_var.rds")
sample0<-readRDS("acral_sub_rppa/ds/sample_info0.rds")




sig_abs<-stats0%>%
  filter(t_test_p_excl<=0.05)%>%
  pull(antibody_name)



mat1<-mat0[sig_abs,]


zmat1<-t(scale(t(mat1)))


hm1<-Heatmap(zmat1)

