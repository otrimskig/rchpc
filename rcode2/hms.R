setwd("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig")

library(tidyverse)

library(ComplexHeatmap)

library(EnhancedVolcano)









# de <- read_csv("de-nf1 KO; pten KO; ink KO; atrx wt vs. nf1 KO; pten KO; ink KO; atrx KO.csv")
# 
# hm<-de%>%
#   filter(!is.na(gene_name))%>%
#   
#   filter(abs(logFC)>=1.5)%>%
#   
#   
#   filter(FDR<.1)%>%
#   
#   group_by(gene_name)%>%
#   slice(1)%>%
#   ungroup()%>%  
#   column_to_rownames("gene_name")%>%
#   
#   select(starts_with("rpkm"))
# 
# 
# 
# 
# hp<-Heatmap(as.matrix(hm))
# 
# sample_info2<-sample_info%>%
#   mutate(col_name=paste0("rpkm_", sample_id))





all<-readRDS("23908R/v05-all_counts_plus_info.rds")%>%
  
  rename(gene_name=external_gene_name)%>%
  filter(!is.na(gene_name))%>%
  
  
  #filter(abs(logFC)>=1.5)%>%
  
  
  #filter(FDR<.1)%>%
  
  group_by(mouse_num, gene_name)%>%
  slice(1)%>%
  ungroup()


all2<-all%>%
  select(gene_name, mouse_num, rpkm)%>%
  mutate(rpkm=log(rpkm+1))%>%
  pivot_wider(values_from = rpkm, names_from = mouse_num)%>%
  
  column_to_rownames("gene_name")%>%
  
  sample_n(1000)
  




all_hm1<-Heatmap(as.matrix(all2), show_row_names = FALSE)



annos<-all%>%
  group_by(mouse_num)%>%
  slice(1)%>%
  ungroup()%>%
  select(mouse_num, resultant_geno, patho_cat, patho_cat_det, aod)%>%
  
  rename(Sample = mouse_num)%>%
  arrange(Sample)%>%
  mutate(index=1:n())%>%
  
  left_join(tibble(index=column_order(all_hm1))%>%
              mutate(order=1:n())
  )%>%
  arrange(order)%>%
  select(resultant_geno, patho_cat, patho_cat_det, aod)


gs<-annos%>%group_by(resultant_geno)%>%slice(1)%>%pull()

all_hm2<-Heatmap(as.matrix(all2), show_row_names = FALSE,
                 top_annotation=HeatmapAnnotation(df=annos, 
                                                  col= list(resultant_geno=c("nf1 KO; pten KO; ink KO; atrx KO"="#66c2a5", 
                                                                             "nf1 KO; pten KO; ink KO; atrx wt"="#fc8d62", 
                                                                             "nf1 wt; pten KO; ink KO; atrx KO"="navy", 
                                                                             "nf1 wt; pten wt; ink wt; atrx wt"="#e7298a")
                                                            ),
                                                  gp = gpar(col = "black")),
                 
                 )

draw(all_hm2)







de1<-read_csv("23908R/de_comps/de-nf1 wt; pten wt; ink wt; atrx wt vs. nf1 KO; pten KO; ink KO; atrx KO.csv")


dem<-de1%>%
  filter(!is.na(gene_name))%>%
  filter(abs(logFC)>2)%>%
  filter(FDR<.0001)%>%
  
  
  arrange(desc(logFC))%>%
  group_by(gene_name)%>%slice(1)%>%ungroup()%>%
  
  select(gene_name, starts_with("rpkm"))%>%
  column_to_rownames("gene_name")


#Heatmap(as.matrix(dem), cluster_rows = FALSE)


library(ggplot2)
library(forcats)


p<-de1%>%
filter(!is.na(gene_name))%>%
  filter(abs(logFC)>2)%>%
  filter(FDR<.000000001)%>%
  
  mutate(gene=fct_reorder(gene_name, logFC))%>%
  mutate(myc=if_else(logFC>0, "y", "n"))

ggplot(p, aes(y=gene, x=0, color=myc))+
  geom_segment(aes(xend=logFC, yend=gene))+
  geom_point(aes(x=logFC),size=1)+
 
  theme_classic()+
  labs(x="logFC")+
  ggtitle("wt vs. 4KO at 1E-9 FDR cutoff")+
theme(legend.position = "none")





atrx<-read_csv("23908R/de_comps/de-nf1 KO; pten KO; ink KO; atrx wt vs. nf1 KO; pten KO; ink KO; atrx KO.csv")

atg<-atrx%>%filter(!is.na(gene_name))%>%
  group_by(gene_name)%>%slice(1)%>%ungroup()%>%
  column_to_rownames("gene_name")
  #filter(abs(logFC)>2)%>%
  #filter(FDR<.01)




EnhancedVolcano(atg, 
                lab=rownames(atg),
                x="logFC",
                y="PValue",
                labSize = 5.0,
                drawConnectors = TRUE,
                boxedLabels = TRUE,
                title="3ko-atrx-wt vs 4KO",
                subtitle = "")





atrx<-read_csv("23908R/de_comps/de-nf1 wt; pten wt; ink wt; atrx wt vs. nf1 KO; pten KO; ink KO; atrx KO.csv")

atg<-atrx%>%filter(!is.na(gene_name))%>%
  group_by(gene_name)%>%slice(1)%>%ungroup()%>%
  column_to_rownames("gene_name")
#filter(abs(logFC)>2)%>%
#filter(FDR<.01)

EnhancedVolcano(atg, 
                lab=rownames(atg),
                x="logFC",
                y="PValue",
                labSize = 5.0,
                drawConnectors = TRUE,
                boxedLabels = TRUE,
                title="wt vs 4KO",
                subtitle = "")




atrx<-read_csv("23908R/de_comps/de-nf1 wt; pten wt; ink wt; atrx wt vs. nf1 KO; pten KO; ink KO; atrx wt.csv")

atg<-atrx%>%filter(!is.na(gene_name))%>%
  group_by(gene_name)%>%slice(1)%>%ungroup()%>%
  
  
  mutate(gene_name=tolower(gene_name))%>%
  filter(gene_name=="pten"|gene_name=="cdkn2a"|gene_name=="nf1"|gene_name=="atrx")%>%
  
  column_to_rownames("gene_name")
#filter(abs(logFC)>2)%>%
#filter(FDR<.01)

EnhancedVolcano(atg, 
                lab=rownames(atg),
                x="logFC",
                y="PValue",
                labSize = 5.0,
                drawConnectors = TRUE,
                boxedLabels = TRUE,
                title="wt vs 3KO-atrx-wt",
                subtitle = "")














all%>%
  #select(gene_name, mouse_num, rpkm)%>%
  mutate(gene_name=tolower(gene_name))%>%
  filter(gene_name=="pten"|gene_name=="cdkn2a"|gene_name=="nf1"|gene_name=="atrx"|gene_name=="cdkn2b")%>%
  group_by(resultant_geno)%>%
  
  ggplot(aes(x=resultant_geno, y=rpkm, color=resultant_geno))+
  geom_point(stat="identity", position = "dodge", size=4)+
  #geom_point(stat="identity", position = "dodge", size=4)+
  facet_wrap(~gene_name)+
  theme_classic()


