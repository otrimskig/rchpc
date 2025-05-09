source("libs.R")

data<-readRDS("acral_paired/dexps/dexp-sample_type-acral v. subq.rds")

samples<-readRDS("acral_paired/ds/v00-sample_info.rds")



# name_map<-data%>%
#   select(starts_with("rpkm"))%>%colnames(.)%>%
#   tibble("name1"=.)

         

name_map<-samples%>%
  select(mouse_num, sample_type, sample_id)%>%
  mutate(name1=paste0("rpkm_", sample_id))%>%
  mutate(name2=paste0("x", mouse_num, "_", sample_type))%>%
  select(name1, name2)




data_renamed <- data %>%
  rename(!!!setNames(name_map$name1, name_map$name2))



saveRDS(data_renamed, "acral_paired/outs/mr-acral_paired_rpkms.rds")


writexl::write_xlsx(data_renamed, "acral_paired/outs/mr-acral_paired_rpkms.xlsx")









all_counts<-readRDS("acral_paired/ds/v00-all_counts.rds")


gene_stats<-readRDS("acral_paired/ds/gene_stats.rds")




name_map2<-samples%>%
  select(mouse_num, sample_type, sample_id)%>%
  mutate(name1=paste0("rpkm_", sample_id))%>%
  mutate(name2=paste0("x", mouse_num, "_", sample_type))%>%
  select(-name1)%>%
  rename(name1=sample_id)%>%
  mutate(name1=toupper(name1))%>%
  mutate(name1=sub("X", "x", name1))
  




raw0<-all_counts%>%
  rename(gene_id_ms=GeneID, gene_length=Length)%>%
  left_join(gene_stats)%>%
  relocate(gene_id_ms, gene_name_ms, gene_name_hu, gene_length,
           mean_rpkm, max_rpkm, sd_rpkm)%>%
  rename(!!!setNames(name_map2$name1, name_map2$name2))





write_csv(raw0, "acral_paired/outs/mr-acral_paired_raw_counts.csv")
