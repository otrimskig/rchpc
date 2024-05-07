library(tidyverse)
library(edgeR)



df<-readRDS("23908R/v03-rpkms.rds")

sample_info<-readRDS("23908R/sample_mouse_num.rds")









modelm<-sample_info%>%
  group_by(resultant_geno)%>%
  ungroup()%>%
  filter(resultant_geno == "nf1 KO; pten KO; ink KO; atrx KO"| resultant_geno == "nf1 wt; pten wt; ink wt; atrx wt")%>%
  
  mutate(resultant_geno=as.factor(resultant_geno))%>%
  
  
  pull(resultant_geno)



design2<-model.matrix(~modelm)%>%as_tibble()



counts2<-df%>%
  select(gene_id, sample_name, raw_read_count)%>%
  left_join(sample_info%>%select(c("sample_name"="sample_id", resultant_geno)))%>%
  
  

filter(resultant_geno == "nf1 KO; pten KO; ink KO; atrx KO"| resultant_geno == "nf1 wt; pten wt; ink wt; atrx wt")%>%
select(-resultant_geno)%>%
  
  pivot_wider(names_from = "sample_name", values_from = "raw_read_count")













counts2.c<-counts2%>%select(-gene_id)

counts2.g<-counts2%>%select(gene_id)



dge2<-DGEList(counts=as.matrix(counts2.c), group=design2$`modelmnf1 wt; pten wt; ink wt; atrx wt`,  genes = counts2.g)




keep <- filterByExpr(dge2)

#filter out some genes. 
dge2 <- dge2[keep, , keep.lib.sizes=FALSE]


#calculate normalization for gene counts. 
dge2<-calcNormFactors(dge2)




dge2<-estimateDisp(dge2, design = design2)

exacDGE<-exactTest(dge2, pair = 1:2, dispersion = "auto", rejection.region = "doubletail")

top<-topTags(exacDGE, n=100000)

df2<-as.data.frame(top)






genes_table<-df%>%
  select(-total_reads)%>%
  pivot_wider(names_from =sample_name, names_glue = "{sample_name}_{.value}",  values_from=c(raw_read_count, rpkm))
  left_join(df2)%>%
  
  relocate(c(logFC:last_col()), .after="unique_gene_names")%>%
  
  filter(!is.na(logFC))



#write_csv(genes_table, "genes_table.csv")







# stdedger<-function(DGE=DGE, design = design, mergefile = mergefile, byx = byx, byy = byy)
#   
# {
#   
#   require(edgeR)
#   
#   DGE<-estimateDisp(DGE, design = design)
#   
#   exacDGE<-exactTest(DGE, pair = 1:2, dispersion = "auto", rejection.region = "doubletail")
#   
#   top<-topTags(exacDGE, n=100000)
#   
#   df<-as.data.frame(top)
#   
#   
#   
# }




# df<-merge(df, mergefile, by.x=byx, by.y=byy)
# 
# df




rpkms_sample_info<-df%>%
  rename(sample_id=sample_name)%>%
  left_join(sample_info)

# saveRDS(rpkms_sample_info, "23908R/v04-rpkms_sample_info.rds")






