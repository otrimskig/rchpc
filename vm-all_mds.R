library(tidyverse)
library(ggplot2)
library(ggrepel)
#library(viridis)
library(edgeR)


vm02<-readRDS("ds/vm-02-filtered_rpkms.rds")

sample_info<-readRDS("ds/v07-per_sample_info.rds")


mat<-vm02%>%
  select(gene_id_ms, sample_id, rpkm)%>%
  pivot_wider(names_from = sample_id, values_from = rpkm)%>%
  column_to_rownames("gene_id_ms")%>%
  as.matrix.data.frame()

lib_size<-vm02%>%
  select(sample_id, read_count)%>%
  ungroup()%>%
  group_by(sample_id)%>%
  summarise(lib_size=sum(read_count))%>%
  ungroup()%>%
  arrange(sample_id)%>%
  pull(lib_size)

dge<-DGEList(counts = mat, 
        lib.size = lib_size)



p<-plotMDS.DGEList(dge)


saveRDS(p, "ds/mds_dge.rds")




stop("done")
dfp<-mds1%>%left_join(all_info)





ggplot(dfp, aes(x, y, color=aod, fill=aod))+
  geom_point(shape=16, alpha=.5, size=4)+
  geom_point(shape=1, size=4)+
   
  
  ggtitle("age of death")+
  viridis::scale_color_viridis()+
  viridis::scale_fill_viridis()+
  theme_classic()




ggplot(dfp, aes(x, y, color=tumor_grouping1, fill=tumor_grouping1))+
  geom_point(shape=16, alpha=.5, size=4)+
  geom_point(shape=1, size=4)+
  
  
  ggtitle("first tumor grouping")+
  #viridis::scale_color_viridis()+
  #viridis::scale_fill_viridis()+
  theme_classic()




ggplot(dfp, aes(x, y, color=patho_cat, fill=patho_cat))+
  geom_point(shape=16, alpha=.5, size=4)+
  geom_point(shape=1, size=4)+
  
  
  ggtitle("major pathologist categorization")+
  #viridis::scale_color_viridis()+
  #viridis::scale_fill_viridis()+
  theme_classic()







ggplot(dfp, aes(x, y, color=patho_cat_det, fill=patho_cat_det))+
  
  
  geom_text_repel(aes(label = ifelse(abs(x)>2.7,mouse_num,"")),
                  min.segment.length = 0.5,
                  
                  segment.color = "grey80",
                  force=30,
                  point.padding=1,
                  label.padding=4) +
  
  geom_point(size=4,shape=16, alpha=.5)+
  geom_point(size=4, shape=1)+
  
  
  ggtitle("pathologist subcategorization")+
  #viridis::scale_color_viridis()+
  #viridis::scale_fill_viridis()+
  
  #geom_text(aes(label=ifelse(abs(x)>2.7,as.character(mouse_num),'')))+
  

  
  theme_classic()












ggplot(dfp, aes(x, y, color=patho_cat_det, fill=patho_cat_det))+
  geom_text_repel(aes(label = ifelse(abs(x)>2.7,mouse_num,"")),
                  min.segment.length = 0.5,
                  
                  segment.color = "grey80",
                  force=30,
                  point.padding=2,
                  label.padding=30) +
  
  
  
  geom_point(aes(size=(150-aod+1)), shape=16, alpha=.5)+
  geom_point(aes(size=(150-aod+1)), shape=1)+
  
  scale_size(range = c(1,10))+
  ggtitle("pathologist subcategorization")+
  #viridis::scale_color_viridis()+
  #viridis::scale_fill_viridis()+
  theme_classic()









ggplot(dfp, aes(x, y, color=resultant_geno, fill=resultant_geno))+
  geom_point(aes(size=(150-aod+1)), shape=16, alpha=.5)+
  geom_point(aes(size=(150-aod+1)), shape=1)+
  
  scale_size(range = c(1,10))+
  ggtitle("resultant genotype")+
  #viridis::scale_color_viridis()+
  #viridis::scale_fill_viridis()+
  theme_classic()






ggplot(dfp, aes(x, y, color=resultant_geno, fill=resultant_geno))+
  geom_point(shape=16, alpha=.5, size=4)+
  geom_point(shape=1, size=4)+
  
  
  ggtitle("resultant genotype")+
  #viridis::scale_color_viridis()+
  #viridis::scale_fill_viridis()+
  theme_classic()
