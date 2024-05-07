setwd("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig")
library(tidyverse)

colors<-readRDS("colors.rds")


variable_names<-colors%>%
  filter(!grepl("^aod", variable))%>%
  group_by(variable)%>%slice(1)%>%ungroup()%>%select(variable)%>%pull()

hm_colors_list<-list()
for (i in 1:length(variable_names)){


v<-colors%>%filter(variable==variable_names[i])

hm_colors_list[[variable_names[i]]]<-setNames(v$hex, v$condition)

}

saveRDS(hm_colors_list, "hm_colors_list.rds")

