library(tidyverse)
library(dtplyr)
library(broom)


#load previously made data of rpkms, plus all genes.
df<-readRDS("k19mf/ds/vm-02-filtered_rpkms.rds")
gs<-readRDS("k19mf/ds/gene_stats.rds")


library(googlesheets4)
sheet_id<-"https://docs.google.com/spreadsheets/d/1QOX3OAygT4wLH7ZTzX0WVm31cIGbCn9SqYSTfWJpVv8/edit?gid=0#gid=0"
name_of_sheet<-"Sheet1"

############
#authorize user.
gs4_auth(email = "gotrimski@gmail.com")

#read input sheet
sheets_df<-read_sheet(sheet_id, 
                      sheet = name_of_sheet)%>%
  mutate(across(1:last_col(), as.character))%>%
  
  #then replace all "NULL" with NA. 
  mutate(across(1:last_col(), function(x){na_if(x, "NULL")}))%>%
  
  janitor::clean_names()






genes_to_include1<-gs%>%
  semi_join(sheets_df)%>%
  pull(gene_name_ms)



# 
# #set vector of pattern of interests (mouse gene names start with...)
interests<-c("^Il",
             "^Itg",
             "^Ifn",
             "^Cxcl",
             "^Ccl",
             "^Ccr",
             "^Irf")


# #get vector of gene names based on criteria.
genes_to_include2<-gs%>%
  filter(grepl(paste(interests, collapse = "|"), gene_name_ms))%>%
  pull(gene_name_ms)



genes_to_include<-c(genes_to_include1, genes_to_include2)%>%unique()%>%sort()


all_relevent_genes<-tibble("gene_name_ms"=genes_to_include)%>%left_join(sheets_df)%>%
  mutate(in_kircher_ds=if_else(!is.na(group), "1", "0"))%>%
  
  
  mutate(group=if_else(is.na(group)&grepl("^Il", gene_name_ms), "Interleukins", group))%>%
  
  
  mutate(group=if_else(is.na(group)&grepl("^Ccl", gene_name_ms), "CCL Chemokines", group))%>%
  
  
  
  mutate(group=if_else(is.na(group)&grepl("^Cxcl", gene_name_ms), "CXCR Chemokines and Receptors", group))%>%
  
  mutate(group=if_else(grepl("^Ccr", gene_name_ms), "CCL Chemokine Receptors", group))%>%
  mutate(group=if_else(grepl("^Itg", gene_name_ms), "Integrins", group))%>%
  
  mutate(group=if_else(grepl("^Irf", gene_name_ms), "Interferon Regulatory Factors", group))%>%
  
  mutate(group=if_else(grepl("^Ifn", gene_name_ms), "Interferons", group))






saveRDS(all_relevent_genes, "k19mf/ds/immune2-genes_groups.rds")




include_genes<-all_relevent_genes%>%pull(gene_name_ms)



#pull out all rpkms from dataset for those genes.
immune_reads<-df%>%filter(gene_name_ms %in% include_genes)

#save externally.
saveRDS(immune_reads,"k19mf/ds/immune2_rpkms.rds")





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
    mean_rpkm_acral = mean(rpkm[tumor_type == "acral"], na.rm = TRUE),
    mean_rpkm_subq = mean(rpkm[tumor_type == "subq"], na.rm = TRUE),
    
    # Calculate log2 fold change
    log2_fold_change = log2(mean_rpkm_acral / mean_rpkm_subq),
    pseudo_log2_fold_change = log2(mean_rpkm_acral / mean_rpkm_subq+1E-10),
    
    ######### value ~ group ####
    t_test = list(t.test(rpkm ~ tumor_type))) %>%
  
  mutate(p_unadjusted = map_dbl(t_test, ~ .x$p.value)) %>%
  mutate(p_adj_bonferroni = p.adjust(p_unadjusted, method = "bonferroni"),  # Bonferroni correction
         p_adj_fdr = p.adjust(p_unadjusted, method = "fdr"))%>%
  
  select(-t_test)%>%
  
  
  
  
  mutate(fdr_stars=gtools::stars.pval(p_adj_fdr))%>%
  mutate(fdr_stars=trimws(gsub("\\.", "", fdr_stars)))%>%
  mutate(fdr_stars=gsub("^$", "ns", fdr_stars))



results3<-results1%>%
  left_join(results2)



saveRDS(results3,
        
        "k19mf/ds/immune2-acral-v-subq-stats.rds")

