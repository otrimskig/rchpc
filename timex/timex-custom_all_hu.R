source("libs.R")

#cudc rna seq timex analysis 
#timex does not care about comparisons, so any group loaded is fine 
#needs RPKM data in a format where the rownames are the gene names (gene symbol) and the rpkms are a matrix 


library(tidyverse)
library(dtplyr)
library(GSVA)



mat<-readRDS("nf1g/ds/vm-h-01-rpkms_wide_human.rds")%>%
  column_to_rownames("gene_name_hu")%>%
  data.matrix()



gmt_files<-tibble(fn=list.files(path="timex/ds_hu", full.names = T))%>%
  filter(grepl("symbols", fn))%>%
  pull()


qusage::read.gmt(gmt_files[1])


# msig_list<-list()
# 
# all_hu<-for(g in 1:length(gmt_files)){
#   
# name_element<-basename(gmt_files[g])
#   
# msig_list[[name_element]][["gene_signatures"]]<-qusage::read.gmt(gmt_files[g])
# 
# }
# 
# 
# 


msig_list<-list()

all_hu<-for(g in 1:length(gmt_files)){
  
  name_element<-basename(gmt_files[g])
  
  msig_list[[name_element]][["gene_signature"]]<-qusage::read.gmt(gmt_files[g])
  
}




# msig_list<-msig_list[1]

for (na in 1:length(names(msig_list))){

element_name<-names(msig_list)[na]


Signature_list <- msig_list[[element_name]][["gene_signature"]]

ssgeaP<-ssgseaParam(mat, Signature_list)
suppressWarnings(gng_ssgsea<- gsva(ssgeaP))


msig_list[[element_name]][["gsva_values"]]<-gng_ssgsea%>%as.data.frame()%>%as.matrix.data.frame()


}



saveRDS(msig_list, "timex/ds/hu-msig_list-all-gsva-values.rds")


# 
# msig_list<-readRDS("timex/ds/hu-msig_list-all-gsva-values.rds")




# bs<-1

inset<-tibble()

for(bs in 1:length(names(msig_list))){

element_name<-names(msig_list)[bs]


addition<-msig_list[[element_name]][["gsva_values"]]%>%
  as.data.frame()%>%
  rownames_to_column("original_pathway_name")%>%
  tibble()%>%
  mutate(gmt_file=element_name)%>%
  relocate(gmt_file)


inset<-bind_rows(inset,addition)

}



inset0<-inset%>%
  group_by(original_pathway_name)%>%
  arrange(gmt_file)%>%
  slice(1)%>%
  ungroup()
  

inset1 <- inset0 %>%
  mutate(updated_pathway_name = paste0(substr(gmt_file, 1, 2), ".", original_pathway_name))%>%
  relocate(updated_pathway_name)%>%
  select(-gmt_file, -original_pathway_name)%>%
  column_to_rownames("updated_pathway_name")%>%
  as.matrix.data.frame()



saveRDS(inset1, "timex/ds/hu-msig_all-gsva-values.rds")

