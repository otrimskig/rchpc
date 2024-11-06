source("libs.R")
library(tidyverse)
library(dtplyr)
library(ggplot2)
library(ggbeeswarm)
library(RColorBrewer)
library(ggdist)



mn<-readRDS("nf1g/ds/mutations_at_nf1_guide_site.rds")%>%
  mutate(sample_id=tolower(paste0("x", sample)))%>%
  rename(total_mut_reads=total_reads)


sa<-readRDS("nf1g/ds/v10-per_sample_updated.rds")%>%
  left_join(mn)



sa_sub<-sa%>%filter(patho_grade %in% c("2", "3", "4"))

plot<-ggplot(sa_sub, aes(x=mutation_perc, color=resultant_geno, y=factor(patho_grade, levels=c("NED", "N", "2", "3", "4"))))+
  
  geom_beeswarm(size=3, alpha=.7, method="square", cex=5)+
  
  facet_grid(patho_cat_name ~ .,  switch = "y")+
  theme_bw()+
  
  scale_y_discrete(position = "right")+
  scale_x_continuous(limits = c(0,100))+
  
  theme(strip.text.y.left = element_text(angle = 0),
       legend.position ="none")+
  
  
  labs(y="Patho Grade",
       x="Mutation Percentage")+
  theme(axis.title.y=element_text(angle=45))+
  scale_colour_brewer(palette = "Dark2")

# ggplot(sa_sub, aes(x=mutation_perc, y=factor(patho_grade, levels=c("NED", "N", "2", "3", "4"))))+
#   geom_point(aes(color=patho_cat_name), size=3, alpha=.7)+
#   facet_grid(patho_cat_name ~ ., scales = "free_y")+
#   theme_classic()


  # 






ggsave("nf1g/plots/dp-nf1-mutation-perc.pdf",
       plot=plot,
       height=10,
       width=30,
       scale = .5,
       dpi=600,
       limitsize = FALSE
)





