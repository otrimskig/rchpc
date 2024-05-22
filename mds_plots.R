library(tidyverse)
library(ggrepel)

#get all sample metadata
all_info<-readRDS("ds/v10-per_sample_updated.rds")


#get coordinates from plot-mds analysis, and dimension length.
mds<-readRDS("ds/mds_dge.rds")
mds_coords<-tibble(x=mds$x,
       y=mds$y,
       sample_id=attributes(mds[["distance.matrix.squared"]])[["dimnames"]][[1]])
x_dim<-round(mds[["var.explained"]][1]*100, 1)
y_dim<-round(mds[["var.explained"]][2]*100, 1)

dim_ratio<-round(y_dim/x_dim, 2)

#join mds coordinates to sample info.
dfp<-mds_coords%>%left_join(all_info)




colors<-readRDS("ds/hm_colors_list.rds")







colors$patho_cat



#plot 1: pathologist: major categorization, scaled to aod.
plot_patho_cat<-ggplot(dfp, aes(x, y, color=patho_cat_name, fill=patho_cat_name))+
  ggtitle("Pathologist Major Categorization")+
  geom_text_repel(aes(label = mouse_num),
                  min.segment.length = 0,
                  segment.color = "grey80",
                  force=20,
                  point.padding=15,
                  show.legend = FALSE) +
  
  geom_point(aes(size=(150-aod+1)), shape=16, alpha=.5)+
  geom_point(aes(size=(150-aod+1)), shape=1)+
  
  scale_size(range = c(2,10))+
  scale_fill_discrete(colors$patho_cat)+
  
  labs(subtitle = "size of point scaled as days remaining to experimental endpoint upon death")+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  guides(fill="none")+
  guides(size="none")+
  
  labs(x=paste0("Dim X (", x_dim, "% of difference)"))+
  labs(y=paste0("Dim Y (", y_dim, "% of difference)"))+
  
  labs(color="Category")+
  coord_fixed(ratio=dim_ratio)+
  theme_classic()

plot_patho_cat













#plot 2: pathologist sub category, scaled to aod.
plot_patho_subcat<-ggplot(dfp, aes(x, y, color=patho_cat_det_name, fill=patho_cat_det))+
  ggtitle("Pathologist Sub-categorization")+
  
  geom_text_repel(aes(label = mouse_num),
                  min.segment.length = 0,
                  segment.color = "grey80",
                  force=20,
                  point.padding=15,
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
  coord_fixed(ratio=dim_ratio)+
  theme_classic()

#plot_patho_subcat


#plot 3: resultant "genotype" scaled to aod.
plot_resultant_geno<-ggplot(dfp, aes(x, y, color=cohort_proper_name, fill=resultant_geno))+
  ggtitle("Resultant Genotype")+
  
  geom_text_repel(aes(label = mouse_num),
                  min.segment.length = 0,
                  segment.color = "grey80",
                  force=20,
                  point.padding=15,
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
  coord_fixed(ratio=dim_ratio)+
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
       width = 40,
       height =7,
       unit="in"
       
)
}

