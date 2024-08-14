library(tidyverse)

k19<-read_csv("../exp_data/kircher19/featurecounts/kircher19-all_sample_fc_stats.csv")%>%
  rename(srr_id=sample_id)%>%
  left_join(read_csv("kircher19/by_sample-w-srr.csv"))%>%
  mutate(exp="k19")%>%
  mutate(sample_id=gnomex_id)%>%
  mutate(sample_id = if_else(nchar(sample_id)==7, 
                      paste0(substr(sample_id, 1,6), "0", substr(sample_id, 7,7)),
         sample_id)
          )




mf<-read_csv("../exp_data/24385R/featurecounts/24385-all_sample_fc_stats.csv")%>%
  mutate(exp="mf5")%>%
  mutate(sample_id2=sample_id)%>%
  mutate(sample_id=substr(sample_id, 1,7))%>%
  mutate(sample_id=paste0(substr(sample_id, 1,6), "0", substr(sample_id, 7,7)))




