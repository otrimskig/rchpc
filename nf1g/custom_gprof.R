source("libs.R")

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(gprofiler2)

data1<-readRDS("nf1g/gresult1test.rds")


data2 <- data1%>%
  mutate(point_index=1:n())%>%
  rename(pathway_list=source)


p1 <- ggplot(data2, aes(x = pathway_list,
                        y = -log10(p_value), 
                        label = point_index)) +
  # generate jittered points
  geom_jitter(position = position_jitter(width = 0.2, seed = 3))+
  coord_cartesian(clip = "off")
  
  
data_coord <- ggplot_build(p1)$data[[1]]%>%
  select(x,y,label)%>%
  rename(point_index=label, x_coord=x, neg_log10_p_value=y)

data3<-data2%>%
  left_join(data_coord)%>%
  mutate(group_x_coord=round(x_coord, digits = 0))%>%
  mutate(dist_from_group_x=x_coord-group_x_coord)%>%
  mutate(label_x_pos=if_else(dist_from_group_x<0, group_x_coord-.3, group_x_coord+.3))%>%
  
  mutate(g_label = if_else(neg_log10_p_value > 5, 
                           stringr::str_wrap(term_name, width=25),
                           ""))
  
  
  

my_pal <- function(range = c(2, 12)) {
  force(range)
  function(x) scales::rescale(x, to = range, from = c(0, 1))
}







p2 <- ggplot(data3, aes(x = x_coord, 
                        y = neg_log10_p_value, 
                        label = g_label)) +
  
  
  
  
  # Add jittered points
  geom_point(aes(size = intersection_size, 
                 fill = pathway_list),
              color = "black",
              alpha = 0.6, 
              shape = 21,      # Ensures both fill & stroke work
              stroke = 0.5) +  # Controls outline thickness
  
  #facet_wrap(vars(pathway_list), nrow = 1)+

  
  
  geom_label_repel(data=subset(data3),
                   aes(label = g_label,
                       point.size = intersection_size),
                   xlim=c(1.4,5),
                   #direction = "both",
                   point.padding = 0,
                   min.segment.length = 0,
                   #size = 3,
                   # segment.ncp = 3,
                   # segment.curvature = -.1,
                   max.overlaps = 10,
                   force = 2,
                   hjust= 0,
                   box.padding = 0)+
 
   # geom_label_repel(data=subset(data3, dist_from_group_x >= 0),
   #                 aes(label = g_label,
   #                     point.size = intersection_size),
   #                 xlim=c(.5,4.5),
   #                 direction = "both",
   #                 point.padding = 0,
   #                 min.segment.length = 0,
   #                 size = 3,
   #                 # segment.ncp = 3,
   #                 # segment.curvature = -.1,
   #                 max.overlaps = 5,
   #                 force = 1,
   #                 hjust= 0,
   #                 box.padding = 0) +
  
  # 

  
  
  
  





  continuous_scale(
    aesthetics = c("size", "point.size"),  # Scale both point size and label size
    palette = my_pal(c(1, 10))  # Custom palette for scaling
  ) +
  coord_cartesian(clip = "off")+
  theme_classic()




print(p2)
