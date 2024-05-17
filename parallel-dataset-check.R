library(tidyverse)
library(dtplyr)


#get vector of all .rds files containing dexp analyses.
all_dexps<-list.files("dexps", full.names = T, pattern="^dexp.*\\.rds$")

results_list <- list()

for (i in 1:length(all_dexps)) {

de<-readRDS(all_dexps[i])%>%
  
  
  
  
  summarize(pv=min(PValue), maxfc=max(logFC), minfc=min(logFC))


results_list[[i]] <- de

}


lp<-purrr::map_dfr(results_list, bind_rows)%>%

summarise(pv=min(pv), maxfc=max(maxfc), minfc=min(minfc))


  
lvp<-c(lp)


saveRDS(lvp, "ds/de_axis_lims.rds")



