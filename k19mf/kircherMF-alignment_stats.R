library(tidyverse)

k19<-read_csv("../exp_data/kircher19/featurecounts/kircher19-all_sample_fc_stats.csv")%>%
  rename(srr_id=sample_id)%>%
  left_join(read_csv("kircher19/by_sample-w-srr.csv"))%>%
  mutate(exp="k19")%>%
  mutate(sample_id=gnomex_id)%>%
  mutate(sample_id = if_else(nchar(sample_id)==7, 
                      paste0(substr(sample_id, 1,6), "0", substr(sample_id, 7,7)),
         sample_id)
          )%>%
  mutate(mouse_num=as.character(mouse_num))




mf<-read_csv("../exp_data/24385R/featurecounts/24385-all_sample_fc_stats.csv")%>%
  mutate(exp="mf5")%>%
  mutate(sample_id2=sample_id)%>%
  mutate(sample_id=substr(sample_id, 1,7))%>%
  mutate(sample_id=paste0(substr(sample_id, 1,6), "0", substr(sample_id, 7,7)))%>%
  relocate(sample_id)




mf_sample<-readxl::read_xlsx("kircher19/Mouse Information-RNA Seq.xlsx")%>%
  janitor::clean_names()%>%
  mutate(sample_id=paste0(substr(gnomex_id, 1,6), "0", substr(gnomex_id, 7,7)))%>%
  relocate(sample_id)


mf2<-mf%>%left_join(mf_sample)%>%
  mutate(mouse_num=as.character(mouse_num))%>%
  rename(genotype = strain, aod=age)




all_ids<-c(k19%>%pull(sample_id), 
       mf2%>%pull(sample_id))%>%
  as_tibble()

all<-bind_rows(k19, mf2)%>%
  relocate(sample_id)%>%
  select(-c(library_protocol:last_col()))%>%
  select(where(~ n_distinct(.) > 1))%>%
  select(-sample_type)%>%
  select(-sample_id2)%>%
  mutate(organism="Mus musculus")



write_csv(all, "kircher19/all_by_sample_stats.csv", na="")

