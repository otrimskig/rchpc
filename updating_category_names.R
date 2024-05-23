library(tidyverse)
library(readxl)
library(dtplyr)

#use filtered rpkms to ensure no information is carried over from now-excluded data.
include<-readRDS("ds/vm-02-filtered_rpkms.rds")%>%group_by(sample_id)%>%slice(1)

sample_info<-readRDS("ds/v07-per_sample_info.rds")%>%
  semi_join(include, by="sample_id")



xl_path<-"ds/nf1-renaming dataset.xlsx"

sheets<-excel_sheets(xl_path)



for(i in 1:length(sheets)){
  
x<-read_excel(xl_path, sheet=sheets[i])

new_info_col<-colnames(x)[2]

factor_v<-x%>%pull(2)

x2<-x%>%
  mutate(!!sym(new_info_col) := factor(!!sym(new_info_col), levels=factor_v))





sample_info<-left_join(sample_info, x2)

}

sample_info2<-sample_info%>%
  select(sort(names(.)))%>%
  relocate(patho_cat_name, .after="patho_cat")%>%
  relocate(patho_cat2_name, .after="patho_cat2")%>%
  relocate(patho_cat_det_name, .after="patho_cat_det")%>%
  relocate(sample_id, mouse_num, total_reads)






saveRDS(sample_info2, "ds/v10-per_sample_updated.rds")




