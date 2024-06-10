library(tidyverse)
library(gprofiler2)




genes<-readRDS("ds/gene_stats.rds")

#all data. 
gostres<- gost(query = genes$gene_id_ms,
               organism = "mmusculus",
               numeric_ns = "ENTREZGENE_ACC")


gostplot(gostres, capped = TRUE, interactive = TRUE)


#genericise for m-dexps. 


dexp_files<-list.files(path = "m-dexps", pattern= ".rds$", full.names = TRUE)


for (i in 1:length(dexp_files)){
  
i<-1

d<-readRDS(dexp_files[i])%>%
  filter(FDR<0.01)%>%
  filter(abs(logFC)>2)%>%
  arrange(logFC)


d_name<-basename(dexp_files[i])



g<- gost(query = d$gene_id_ms,
               organism = "mmusculus",
               numeric_ns = "ENTREZGENE_ACC",
         ordered_query = TRUE)

p<-gostplot(g, capped = TRUE, interactive = F)

}


p



publish_gosttable(g)
