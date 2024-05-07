library(tidyverse)
library(janitor)


setwd("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/23908R")

df<-readRDS("v03-rpkms.rds")


samples<- read_delim("23908R-samples.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)%>%
  select(ID, `Sample Name`)%>%
  mutate(sample_no = sprintf("%02d", as.numeric(sub(".*X(.*)", "\\1", ID))))%>%
  mutate(sample_id = paste0("x23908x", sample_no))%>%
  
  
  mutate(mouse_num = str_extract(`Sample Name`, "\\b\\d{5}\\b"))



info2<-read_csv("sample_info2.csv")%>%
  janitor::clean_names()%>%
  select(mouse_num, tumor, aod, resultant_geno)%>%
  rename(tumor_grouping1=tumor)%>%
  mutate(mouse_num = as.character(mouse_num))



sample2<-samples%>%
  left_join(info2, by="mouse_num")%>%
  select(-c(`ID`, `Sample Name`, sample_no))%>%
  
  
  mutate(resultant_geno = if_else(mouse_num=="36183"|mouse_num=="36184", "nf1 wt; pten wt; ink wt; atrx wt", resultant_geno))

saveRDS(sample2, "sample_mouse_num.rds")



df2<-df%>%
  rename(sample_id=sample_name)%>%
  left_join(sample2, by="sample_id")



