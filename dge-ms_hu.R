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

mousegenes<-dgetable$external_gene_name



dgetable4<-convertMousetoHuman(mousegenes, dgetable)