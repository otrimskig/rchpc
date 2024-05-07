library(tidyverse)
library(edgeR)

setwd("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig")

df<-readRDS("23908R/v04-rpkms_sample_info.rds")




counts<-df%>%
  
  filter(sample_id!="x23908x21"&sample_id!="x23908x15")%>%
  
  select(gene_id, mouse_num, raw_read_count)%>%

  
  
  pivot_wider(names_from = "mouse_num", values_from = "raw_read_count")

counts.c<-counts%>%
  select(-gene_id)

counts.g<-counts%>%
  select(gene_id)


counts.f<-df%>%
  filter(sample_id!="x23908x21"&sample_id!="x23908x15")%>%
  mutate(group=as.factor(resultant_geno))%>%
  group_by(sample_id)%>%
  slice(1)%>%
  ungroup()%>%
  select(group)








dge<-DGEList(counts=counts.c, group=counts.f$group, genes=counts.g)

x<-nlevels(counts.f$group)

gindex<-counts.f%>%
  group_by(group)%>%
  slice(1)%>%
  ungroup()%>%
  mutate(gindex=1:n())


gindex2<-counts.f%>%
  left_join(gindex)






p<-plotMDS(dge, col=gindex2$gindex)


#extracting plottable info from MDS calculation. 
dfp<-tibble(sample_id=colnames(as_tibble(p$distance.matrix.squared)),
            x=p$x, y=p$y)


saveRDS(dfp, "mds.rds")



ggplot(dfp, aes(x, y, color=sample_id))+
  geom_point()

