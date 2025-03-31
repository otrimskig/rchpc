source("libs.R")
library(tidyverse)
library(fs)



gsva_info<-tibble(path=dir_ls("nf1g/ds/gsva"))%>%
  filter(grepl("^nf1g/ds/gsva/gsva-u", path))%>%
  mutate(bn=basename(path))%>%
  mutate(pathway_list=sub("^gsva-u-", "", bn))%>%
  mutate(pathway_list=sub(".rds$", "", pathway_list))

a2<-tibble()

 
for (p in 1:length(gsva_info$path)){



a0<-readRDS(gsva_info$path[p])


a1<-a0%>%
  as.data.frame() %>%
  rownames_to_column(var = "rowname")%>%
  rename(pathway=rowname)%>%
  
  mutate(pathway_list=gsva_info$pathway_list[p])%>%
  relocate(pathway_list)


a2<-bind_rows(a1,a2)



}
















a3 <- a2 %>%
  group_by(pathway)%>%
  arrange(pathway_list)%>%
  slice(1)%>%
  ungroup()%>%
  arrange(pathway_list, pathway) %>%
  mutate(pathway_index =  paste0("m", sprintf("%05d", 1:n()))) %>%
  relocate(pathway_index) %>%
  rename(pathway_list_file = pathway_list)





b1<-a3%>%
  select(1:3)%>%

  mutate(pathway_code_1 = sub("\\..*", "", pathway_list_file))%>%
  
  mutate(pathway_code_x = if_else(grepl("_", pathway), sub("_.*", "", pathway), NA))%>%
  
  
  mutate(pathway_code_1_desc= if_else(pathway_code_1=="m1", "positional gene sets", 
                                      
                                      if_else(pathway_code_1=="m2", "curated gene sets",
                                              
                                          if_else(pathway_code_1=="m3", "regulated gene sets",
                                                  
                                                  if_else(pathway_code_1=="m5", "ontology gene sets",
                                                          
                                                          if_else(pathway_code_1=="m8", "cell type signature gene sets", 
                                                                  
                                                                  if_else(pathway_code_1=="mh", "orthology-mapped hallmark gene sets",

                                                                  NA)))))))%>%
  
  
  mutate(pathway_code_2_desc= if_else(pathway_code_1=="m2", 
                                      if_else(pathway_code_x=="BIOCARTA"|
                                                pathway_code_x=="REACTOME"|
                                                pathway_code_x=="WIKIPATHWAYS", "canonical pathways", "chemical and genetic perturbations"), NA))%>%



  mutate(pathway_code_2= if_else(pathway_code_1=="m2", 
                                    if_else(pathway_code_x=="BIOCARTA"|
                                              pathway_code_x=="REACTOME"|
                                              pathway_code_x=="WIKIPATHWAYS", "CP", "CGP"), NA))%>%
  
  
  mutate(pathway_code_2_desc=if_else(pathway_code_1=="m3",
                                if_else(pathway_code_x=="MIR"|pathway_code_x=="LET", 
                                        "miRDB microRNA targets", "GTRD transcription factor targets"), pathway_code_2_desc))%>%
  
  
  mutate(pathway_code_2=if_else(pathway_code_1=="m3",
                                if_else(pathway_code_x=="MIR"|pathway_code_x=="LET", 
                                        "MIRDB", "GTRD"), pathway_code_2))%>%
  
  
  
  mutate(pathway_code_2=if_else(pathway_code_1=="m5", 
                                if_else(grepl("^GO", pathway_code_x), "GO", "MPT"),
                                  pathway_code_2))%>%
  
  
  mutate(pathway_code_2_desc=if_else(pathway_code_1=="m5", 
                                if_else(grepl("^GO", pathway_code_x), "Gene Ontology", "Mouse Phenotype Ontology MP Tumor"),
                                pathway_code_2_desc))%>%
  
  
  mutate(pathway_code_3=if_else(pathway_code_x=="GOBP", "GO:BP", 
                                if_else(pathway_code_x=="GOCC", "GO:CC",
                                        if_else(pathway_code_x=="GOMF", "GO:MF",
                                                if_else(pathway_code_x=="BIOCARTA", "CP:B",
                                                        if_else(pathway_code_x=="REACTOME", "CP:REACTOME",
                                                                if_else(pathway_code_x=="WP", "CP:WIKIPATHWAYS", NA))))))
         )%>%

  mutate(pathway_code_3_desc=if_else(pathway_code_x=="GOBP", "GO biological process", 
                              if_else(pathway_code_x=="GOCC", "GO cellular component",
                                      if_else(pathway_code_x=="GOMF", "GO molecular function",
                                              if_else(pathway_code_x=="BIOCARTA", "BioCarta gene sets",
                                                      if_else(pathway_code_x=="REACTOME", "Reactome gene sets",
                                                              if_else(pathway_code_x=="WP", "WikiPathways gene sets", NA))))))
       )%>%
  
  
  rename(pathway_prefix=pathway_code_x, pathway_name=pathway)%>%
  
  relocate(pathway_index, pathway_list_file, pathway_prefix, pathway_name,
           pathway_code_1, pathway_code_1_desc,
           pathway_code_2, pathway_code_2_desc,
           pathway_code_3, pathway_code_3_desc
           )%>%
  
  
  mutate(pathway_full_name = paste0(pathway_code_1, "-",
                                    pathway_code_2, "-",
                                    pathway_code_3, "-",
                                    pathway_name))%>%
  
  mutate(pathway_full_name= gsub("NA-", "", pathway_full_name))%>%
  
  mutate(pathway_full_name=sub("GO-GO", "GO", pathway_full_name))%>%
  mutate(pathway_full_name=sub("CP-CP", "CP", pathway_full_name))



















#load rda for pathways.
load("timex/ds/allSignatures.rda")  

timex_names<-names(TIMEx)
hallmark_names<-names(Hallmark)%>%gsub("HALLMARK_", "", .)
immune_names<-names(Immune_sig)


tm1<-readRDS("nf1g/ds/gsva_u-onco.rds")



tm2<-rownames(tm1)%>%
  tibble("original_pathway_name"=.)%>%
  mutate(parent_pathway=if_else(grepl("^KEGG", original_pathway_name), "KEGG", NA))%>%
  mutate(parent_pathway=if_else(grepl("^onco_", original_pathway_name), "onco", parent_pathway))%>%
  mutate(parent_pathway = if_else(is.na(parent_pathway) & original_pathway_name %in% timex_names, "TIMEx", parent_pathway))%>%
  mutate(parent_pathway = if_else(is.na(parent_pathway) & original_pathway_name %in% immune_names, "Immune_Sig", parent_pathway))%>%
  mutate(parent_pathway = if_else(is.na(parent_pathway) & original_pathway_name %in% hallmark_names, "Hallmark", parent_pathway))%>%
  
  mutate(naked_names=gsub("^KEGG_", "", gsub("^onco_", "", original_pathway_name)))%>%
  
  mutate(pathway_full_name=paste0("hu-", parent_pathway, "-", naked_names))%>%
  unique()%>%
  mutate(pathway_index = paste0("h", sprintf("%03d", 1:n())))%>%

  rename(pathway_name=original_pathway_name)%>%
  select(-naked_names)%>%
  rename(pathway_code_1_desc=parent_pathway)%>%
  mutate(pathway_list_file="allSignatures.rda")







tm3<-bind_rows(b1, tm2)

samples<-readRDS("nf1g/ds/v10-per_sample_updated.rds")%>%
  select(sample_id,mouse_num)%>%
  mutate(mouse_num=paste0("x", mouse_num))


tm4<-tm1%>%as_tibble(rownames = "pathway_name")%>%
  janitor::clean_names()%>%
  select(pathway_name, samples$mouse_num)%>%
  rename_with(~ samples$sample_id[match(.x, samples$mouse_num)], .cols = samples$mouse_num)%>%
  left_join(tm2)%>%
  select(starts_with("x"), pathway_full_name)%>%
  relocate(pathway_full_name)









all_pathway_meta<-tm3%>%
  relocate(pathway_full_name)



saveRDS(all_pathway_meta, "nf1g/gsvas/ds/gsva_pathways_meta.rds")




b2<-a3%>%
  left_join(b1)%>%
  select(starts_with("x"), pathway_full_name)


all_pathway_gsva<-bind_rows(tm4, b2)%>%
  column_to_rownames("pathway_full_name")%>%
  as.matrix.data.frame()



saveRDS(all_pathway_gsva, "nf1g/gsvas/ds/gsva_pathways_matrix.rds")
