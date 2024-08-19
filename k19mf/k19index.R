library(tidyverse)
library(dtplyr)

k1<-data.table::fread("k19mf/ds/GSE122781_star_counts22.txt")

k2<-data.table::fread("k19mf/ds/GSE122781_series_matrix.txt", fill=TRUE)

library(janitor)


meta0<-k2%>%
  slice(1:31)%>%
  t()%>%
  row_to_names(1)%>%clean_names()%>%
  as_tibble()%>%
  slice(1)

rownames(meta0) <-NULL


write_csv(meta0, "k19mf/ds/meta.csv")


###sample data
sample0<-k2%>%
  slice(32:n())%>%
  t()%>%as_tibble()%>%
  row_to_names(1)%>%clean_names()

sample1<-sample0%>%
  select(-sample_status)%>%
  select(-3)%>%
  select(-3)%>%
  select(-4)%>%
  select(-starts_with("sample_extr"))%>%
  select(-starts_with("sample_data_pr"))%>%
  select(-starts_with("series_matrix"))%>%
  select(-starts_with("sample_supp"))%>%
  select(-starts_with("sample_cont"))%>%
  select(-starts_with("sample_libr"))%>%
  select(-starts_with("sample_data_row"))%>%
  select(-starts_with("sample_taxid"))%>%
  select(-starts_with("sample_molecule"))%>%
  select(-starts_with("sample_type"))%>%
  select(-"sample_characteristics_ch1_2")%>%
  
  #rename(geo_accession = sample_geo_accession)%>%
  rename(source_type = sample_source_name_ch1)%>%
  
  rename(organism = sample_organism_ch1)%>%
  
  rename(gnomex_id=sample_description_2)%>%
  relocate(gnomex_id)%>%

  rename(mouse_num = sample_title)%>%
  relocate(mouse_num)%>%
  
  
  rename(strain = sample_characteristics_ch1)%>%
  mutate(strain = "Dct::TVA")%>%
  
  rename(aod = sample_characteristics_ch1_3)%>%
  mutate(aod = as.numeric(substr(aod, nchar(aod)-1, nchar(aod))))%>%
  
  
  rename(genotype=sample_characteristics_ch1_4)%>%
  mutate(genotype = sub("genotype: ", "", genotype))%>%
  
  
  mutate(exclusion = if_else(grepl("QC issue", sample_description), "1", NA_character_))%>%
  select(-sample_description)%>%
  
  mutate(sra_id = sub(".*=", "", sample_relation_2))%>%
  relocate(sra_id, .after = "gnomex_id")%>%
  select(-sample_relation_2)%>%
  
  
  mutate(biosample = substr(sample_relation, 51, nchar(sample_relation)))%>%
  select(-sample_relation)%>%
  
  mutate(genotype= paste0(strain, "; ", genotype))%>%
  select(-strain)%>%
  
  relocate(genotype, aod, .after = "mouse_num")


write_csv(sample1, "k19mf/ds/by_sample.csv")
