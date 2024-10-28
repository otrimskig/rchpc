library(tidyverse)
library(dtplyr)
library(broom)


#load previously made data of rpkms, plus all genes.
df<-readRDS("k19mf/ds/vm-02-filtered_rpkms.rds")
gs<-readRDS("k19mf/ds/gene_stats.rds")



#set vector of pattern of interests (mouse gene names start with...)
interests<-c("^Il",
"^Ifn",
"^Cxcl",
"^Ccl",
"^Cd274",
"^Pdcd",
"^Itg")


#get vector of gene names based on criteria.
genes_to_include<-gs%>%
  filter(grepl(paste(interests, collapse = "|"), gene_name_ms))%>%
  pull(gene_name_ms)


#pull out all rpkms from dataset for those genes.
immune_reads<-df%>%filter(gene_name_ms %in% genes_to_include)

#save externally.
saveRDS(immune_reads,"k19mf/ds/immune_reads.rds")





#now calculate stats for all genes. 

#re-name samples acral for future uses. Don't need to run again.

# readRDS("k19mf/ds/vm-00-sample_info.rds")%>%
#   mutate(tumor_type=gsub("foot pad", "acral", tumor_type))%>%
#   saveRDS("k19mf/ds/vm-00-sample_info.rds")





tumor_type<-readRDS("k19mf/ds/vm-00-sample_info.rds")%>%
  select(sample_id, mouse_num, tumor_type)




immune_data_plus_tumor_type<-immune_reads%>%select(gene_name_ms, sample_id, rpkm)%>%
  left_join(tumor_type)%>%
  mutate(tumor_type=gsub("foot pad", "acral", tumor_type))




results1 <- immune_data_plus_tumor_type %>%
  group_by(gene_name_ms,tumor_type) %>%

  summarise(n=n(),
            
            mean=mean(rpkm),
            sd=sd(rpkm),
            
            se=sd/sqrt(n))




results2<-immune_data_plus_tumor_type%>%
  
  group_by(gene_name_ms)%>%
    
    summarise(
    
      ######### value ~ group ####
      t_test = list(t.test(rpkm ~ tumor_type))) %>%
  
  mutate(p_unadjusted = map_dbl(t_test, ~ .x$p.value)) %>%
  mutate(p_adj_bonferroni = p.adjust(p_unadjusted, method = "bonferroni"),  # Bonferroni correction
         p_adj_fdr = p.adjust(p_unadjusted, method = "fdr"))%>%
  
  select(-t_test)
  







results3<-results1%>%
  left_join(results2)



saveRDS(results3,
  
  "k19mf/ds/immune-acral-v-subq-stats.rds")
