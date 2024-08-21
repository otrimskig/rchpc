source("libs.R")

source("k19mf/01.1-bam_to_feature_counts.R")
source("k19mf/01.2-feature_counts-stats.R")
source("k19mf/01.3-counts_to_rds.R")
source("k19mf/01.4-sample_names.R")
source("k19mf/01.5-kircherMF-alignment_stats.R")



library(tidyverse)


counts_check<-readRDS("k19mf/ds/vm-00-all_counts_sample_names_fixed.rds")





counts<-readRDS("k19mf/ds/vm-00-all_counts_sample_names_fixed.rds")%>%
  select(-1, -2)%>%
  summarise_all(sum, na.rm=TRUE)%>%t()%>%
  data.frame()%>%janitor::clean_names()%>%
  
  rownames_to_column()%>%
  rename(sample_id=rowname, total_reads=x)



view(counts_check)
view(counts)



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
  left_join(counts)%>%
  
  
  
  arrange(sample_id)%>%
  mutate(tumor_type=as_factor(tumor_type))


saveRDS(all_info,"k19mf/ds/vm-00-sample_info.rds")




