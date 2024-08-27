source("libs.R")

library(tidyverse)
library(dtplyr)
library(GSVA)



mat<-readRDS("nf1g/ds/ms_ens_rpkms_wide.rds")

list_of_signature_lists<-list.files(path = "timex/ds_ms", full.names = TRUE)


for (i in 1:length(list_of_signature_lists)){


signatures<-qusage::read.gmt(list_of_signature_lists[i])



#make parameters list.
ssGSVA_par<-ssgseaParam(mat, signatures)

#calculate scores
#suppress warnings to ignore errors on sets that are only 1 gene. 
suppressWarnings(ssGSVA_scores<- gsva(ssGSVA_par))



#create unscaled scored object
ssGSVA_u <- t(t(ssGSVA_scores))

saveRDS(ssGSVA_u, 
        paste0("nf1g/ds/gsva/gsva-u-", sub("\\.gmt$", "", basename(list_of_signature_lists[i])), ".rds"))


#create scaled scored object
ssGSVA_z <- t(scale(t(ssGSVA_scores)))

saveRDS(ssGSVA_z, 
        paste0("nf1g/ds/gsva/gsva-z-", sub("\\.gmt$", "", basename(list_of_signature_lists[i])), ".rds"))




}






tb<-tibble(f=list.files("nf1g/ds/gsva", full.names = T))%>%
  mutate(b=basename(f))%>%
  
  mutate(uz=substr(b, 6,6))%>%
  mutate(s=substr(b, 8, nchar(b)))%>%
  mutate(s=sub("\\.rds$", "", s))%>%
  arrange(s)



nested_list <- tb %>%
  group_by(s, uz) %>%
  summarise(f_values = list(f), .groups = 'drop') %>%
  pivot_wider(names_from = uz, values_from = f_values) %>%
  split(.$s) %>%
  map(~ .x %>% select(-s) %>% as.list())


saveRDS(nested_list, "nf1g/ds/gsva_analysis_ds_list_ms.rds")

