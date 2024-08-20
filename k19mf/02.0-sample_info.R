source("libs.R")

library(tidyverse)



#get all sample metadata

###re-do this.
# prev_info<-readRDS("nf1g/ds/v10-per_sample_updated.rds")%>%head()


all_info<-read_csv("k19mf/ds/all_by_sample_stats.csv")%>%
  mutate(genotype=if_else(genotype=="DCT-TVA::Brafca/ca;Ptenf/f;Cdkn2af/f", 
                          "Dct::TVA; BRAFV600E;Cdkn2a-/-;Pten-/-",
                          genotype))%>%
  mutate(sample_id=paste0("x", tolower(sample_id)))%>%
  filter(genotype!="Dct::TVA; BRAFV600E;Cdkn2a-/-")%>%
  mutate(tumor_type=if_else(!is.na(sra_id), "subq", "foot pad"))%>%
  relocate(sample_id, mouse_num, genotype, tumor_type, aod, sex)%>%
  mutate(mouse_num=as.character(mouse_num))%>%
  arrange(sample_id)


saveRDS(all_info,"k19mf/ds/vm-00-sample_info.rds")




