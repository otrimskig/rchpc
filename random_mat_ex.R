

library(tidyverse)
library(dtplyr)
library(GSVA)



mat<-readRDS("rebecca/rpkms_set.rds")%>%
  column_to_rownames("gene_name_hu")%>%
  data.matrix()




fill.matrix = function(expr, nrow=1, ncol=1) {
  matrix(eval(expr, envir=list(x=nrow*ncol)), nrow=nrow, ncol=ncol)
}


x<-mat

mat2<-fill.matrix(rexp(x, rate=.1), nrow=15666, ncol=11)



rownames(mat2)<-rownames(mat)
colnames(mat2)<-v

v<-c()
for (n in 1:length(seq(1:11))){
  v[n]<-c(paste0("x", seq(1:11)[n]))
}



saveRDS(mat2, "rebecca/rpkms_example.rds")
