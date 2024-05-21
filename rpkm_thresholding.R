library(tidyverse)
library(dtplyr)


remove<-as.character(c("x23908x07","x23908x15","x23908x21"))



df<-readRDS("ds/vm-01-gene_id_hu_m_rpkms.rds")%>%
  mutate(max_rpkm = sum(rpkm))%>%
  filter(max_rpkm>=1)%>%
  mutate(sample_id=tolower(sample_id))

df2<-df%>%
  filter(!sample_id %in% remove)%>%
  select(-max_rpkm)


saveRDS(df2, "ds/vm-02-filtered_rpkms.rds")
