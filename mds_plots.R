library(tidyverse)
library(ggrepel)

mds<-readRDS("ds/mds_dge.rds")
all_info<-readRDS("ds/v07-per_sample_info.rds")


mds_coords<-tibble(x=mds$x,
       y=mds$y,
       sample_id=attributes(mds[["distance.matrix.squared"]])[["dimnames"]][[1]])

patho_cat_names<-tibble(patho_cat=c(as.character(c(1,2,3,4)), "NED"),
      patho_cat_name=factor(c("Nerve Sheath Tumors", "Spindle and Epitheliod Tumors", "Gliomas", "Glioneuronal Tumors", "NED"),
                            levels=c("Nerve Sheath Tumors", "Spindle and Epitheliod Tumors", "Gliomas", "Glioneuronal Tumors", "NED")))




dfp<-mds_coords%>%left_join(all_info)%>%
  left_join(patho_cat_names)%>%
  mutate(point_size = 150 - aod + 1,
         point_padding = scales::rescale(point_size, to = c(2, 10)) / 2)
 

x_dim<-round(mds[["var.explained"]][1]*100, 1)
y_dim<-round(mds[["var.explained"]][2]*100, 1)

str(dfp)

#plot 1: pathologist: major categorization, scaled to aod.
plot_patho_cat<-ggplot(dfp, aes(x, y, color=patho_cat_name, fill=patho_cat_name))+
  
  ggtitle("Pathologist Major Categorization")+
  
  geom_text_repel(aes(label = mouse_num),
                 min.segment.length = 0,
                 segment.color = "grey80",
                 force=18,
                 point.padding=12,
                 show.legend = FALSE) +
  
  geom_point(aes(size=(150-aod+1)), shape=16, alpha=.5)+
  geom_point(aes(size=(150-aod+1)), shape=1)+
  
  scale_size(range = c(2,10))+
  
  
  labs(subtitle = "size of point scaled as days remaining to experimental endpoint upon death")+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  guides(fill="none")+
  guides(size="none")+
  
  labs(x=paste0("Dim X (", x_dim, "% of difference)"))+
  labs(y=paste0("Dim Y (", y_dim, "% of difference)"))+
  
  labs(color="Category")+
  theme_classic()


plot_patho_cat

#plot 2: pathologist sub category, scaled to aod.
plot_patho_subcat<-ggplot(dfp, aes(x, y, color=patho_cat_det, fill=patho_cat_det))+
  ggtitle("Pathologist Sub-categorization")+
  
  geom_text_repel(aes(label = mouse_num),
                  min.segment.length = 0,
                  segment.color = "grey80",
                  force=18,
                  point.padding=12,
                  show.legend = FALSE) +
  
  geom_point(aes(size=(150-aod+1)), shape=16, alpha=.5)+
  geom_point(aes(size=(150-aod+1)), shape=1)+
  
  scale_size(range = c(2,10))+
  
  
  labs(subtitle = "size of point scaled as days remaining to experimental endpoint upon death")+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  guides(fill="none")+
  guides(size="none")+
  
  labs(x=paste0("Dim X (", x_dim, "% of difference)"))+
  labs(y=paste0("Dim Y (", y_dim, "% of difference)"))+
  
  labs(color="Category")+
  theme_classic()


#plot 3: resultant "genotype" scaled to aod.
plot_resultant_geno<-ggplot(dfp, aes(x, y, color=resultant_geno, fill=resultant_geno))+
  ggtitle("Resultant Genotype")+
  
  geom_text_repel(aes(label = mouse_num),
                  min.segment.length = 0,
                  segment.color = "grey80",
                  force=18,
                  point.padding=12,
                  show.legend = FALSE) +
  
  geom_point(aes(size=(150-aod+1)), shape=16, alpha=.5)+
  geom_point(aes(size=(150-aod+1)), shape=1)+
  
  scale_size(range = c(2,10))+
  
  
  labs(subtitle = "size of point scaled as days remaining to experimental endpoint upon death")+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  guides(fill="none")+
  guides(size="none")+
  
  labs(x=paste0("Dim X (", x_dim, "% of difference)"))+
  labs(y=paste0("Dim Y (", y_dim, "% of difference)"))+
  
  labs(color="Category")+
  theme_classic()



####################################


to_export<-grep("^plot", ls(), value = TRUE)

for (i in 1:length(to_export)){

name<-to_export[i]
  
object<-get(sym(name))

ggsave(paste0("plots/",fs::path_sanitize(paste0("mds-", name, "-m.pdf"))),
       plot=object,
       
       scale = .75,
       dpi=600,
       width = 15,
       height =10,
       unit="in"
       
)
}

