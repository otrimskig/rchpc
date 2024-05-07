
library(tidyverse)

a<-read_csv("de-nf1 wt; pten wt; ink wt; atrx wt vs. nf1 KO; pten KO; ink KO; atrx KO.csv")%>%
  rename(external_gene_name=gene_name)

#needs mouse genes (dgetable$external_gene_name) and the dgetable

convertMousetoHuman <- function(mousegenes, dgetable){
  
  require("biomaRt")
  
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "dec2021.archive.ensembl.org")
  
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "dec2021.archive.ensembl.org")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = mousegenes , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  dgetable2<-merge(genesV2, dgetable, by.x = "MGI.symbol", by.y ="external_gene_name", all=F, no.dups=T)
  
  
  
  # Print the first 6 genes found to the screen
  
  print(head(dgetable2[,1:3]))
  
  return(dgetable2)
  
}

mousegenes<-a$external_gene_name



dgetable4<-convertMousetoHuman(mousegenes, a)


dt<-dgetable4%>%
  filter(FDR<.1)%>%
  
  select("HGNC.symbol", logFC)%>%
  
  
  
  rename(GENES=1)
  



dt%>%
  as_tibble()%>%
  mutate(logFC=as.numeric(logFC))%>%
  write.table("wt-vs-4KO.txt", row.names = F, col.names = F, quote=F, sep="\t")
  


