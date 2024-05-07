library(tidyverse)
library(edgeR)

setwd("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig")

all_info<-readRDS("23908R/v05-all_counts_plus_info.rds")%>%
  
  #rm x15 and x21 from dataset, the 2 samples that performed poorly in QC. 
  filter(!str_ends(sample_id, "x15|x21"))


all_info%>%
  summarize(sample_id, n())%>%
  view()


gene_names<-all_info%>%
  group_by(gene_id)%>%
  slice(1)%>%
  ungroup()%>%
  select(gene_id, external_gene_name)




models<-all_info%>%
  group_by(sample_id)%>%
  slice(1)%>%
  select(sample_id, mouse_num:last_col())



modelsg<-models%>%group_by(resultant_geno)%>%
  select(resultant_geno)%>%
  slice(1)%>%
  ungroup()%>%
  arrange(desc(resultant_geno))%>%
  mutate(g_index=1:n())


modelsp<-models%>%group_by(patho_cat)%>%
  select(patho_cat)%>%
  slice(1)%>%
  ungroup()%>%
  arrange(desc(patho_cat))%>%
  mutate(p_index=1:n())




modelsd<-models%>%group_by(patho_cat_det)%>%
  select(patho_cat_det)%>%
  slice(1)%>%
  ungroup()%>%
  arrange(desc(patho_cat_det))%>%
  mutate(d_index=1:n())


modelst<-models%>%group_by(tumor_grouping1)%>%
  select(tumor_grouping1)%>%
  slice(1)%>%
  ungroup()%>%
  arrange(desc(tumor_grouping1))%>%
  mutate(t_index=1:n())



models_index<-models%>%
  left_join(modelsd)%>%
  left_join(modelsg)%>%
  left_join(modelsp)%>%
  left_join(modelst)





all_info2<-all_info%>%
  left_join(models_index)






#######################################################################






counts<-all_info2%>%
  ##############change 1/5
  filter(g_index == 3 | g_index== 4)%>%
  arrange(g_index)%>%
  select(sample_id, gene_id, raw_read_count)%>%
  
  pivot_wider(values_from=raw_read_count, names_from = sample_id)


group<-models_index%>%
  arrange(g_index)%>%
  
  ################change 2/5
  filter(g_index == 3 | g_index==4)%>%
  pull(g_index)



design<-model.matrix(~group)


dge<-DGEList(counts=counts%>%select(-gene_id), 
             genes = counts%>%select(gene_id),
             group=group)




keep <- filterByExpr(dge)

#filter out some genes. 
dge <- dge[keep, , keep.lib.sizes=FALSE]


#calculate normalization for gene counts. 
dge<-calcNormFactors(dge)




dge<-estimateDisp(dge, design = design)

exacDGE<-exactTest(dge, pair = 1:2, dispersion = "auto", rejection.region = "doubletail")

top<-topTags(exacDGE, n=100000)

df2<-as.data.frame(top)




df3<-df2%>%
  left_join(gene_names)


file_names<-models_index%>%
  ##########################################change me 3/5
  filter(g_index == 3 | g_index== 4)%>%
  arrange(g_index)%>%
  group_by(g_index)%>%
  slice(1)%>%
  ungroup()%>%
  
  
  ##################var to change 4/5
  pull(resultant_geno)




horizontal_samples<-all_info2%>%
  
  
  
  ####################change me 5/5
  filter(g_index== 3 | g_index ==4)%>%
  arrange(g_index)%>%
  select(sample_id, gene_id, raw_read_count, rpkm)%>%
  
  pivot_wider(values_from = c(raw_read_count, rpkm), names_from = sample_id)


df4<-df3%>%
  rename(gene_name=external_gene_name)%>%
  relocate(gene_name)%>%
  left_join(horizontal_samples)


write_csv(df4, paste0("23908R/de_comps/", "de-",  file_names[1], " vs. ", file_names[2], ".csv"))




