library(tidyverse)




d1_path<-read_delim("sample_meta/Holmen, Sheri Lab_20242024-05-02.txt")%>%
  janitor::clean_names()

d1_seq<-read_delim("sample_meta/Holmen, Sheri Lab_20242024-05-02-seq.txt")%>%
  janitor::clean_names()


d2<-d1_path%>%
  filter(is.na(core_sample_alias))%>%
  select(id, input, sample_name)%>%
  relocate(id)%>%
  rename(path_id=id)%>%
 
  mutate(mouse_num_path = str_extract(sample_name, "\\d{5}"))%>%
  select(-sample_name)



d3<-d1_seq%>%
  select(gnomex_id_mol_diag, id, sample_name)%>%
  mutate(mouse_num_seq = str_extract(sample_name, "\\d{5}"))%>%
  select(-sample_name)%>%
  rename(seq_id = id,
         path_id=gnomex_id_mol_diag)


dj<-left_join(d2,d3)


