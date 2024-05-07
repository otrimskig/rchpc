library(tidyverse)
library(ggplot2)
library(ggrepel)



#input data in tidy format. 
dfp<-mtcars%>%as_tibble()




ggplot(dfp, aes(mpg, hp, size=disp, color=gear, fill=gear))+
  
  #geom for labeling points (from ggrepel)
  #conditional filter for labeling only above certain threshold.
  geom_text_repel(aes(size=60,label = ifelse(mpg>27,mpg,"")),
                
                  
                  #determines layout of label
                  min.segment.length = 0.5,
                  
                  segment.color = "grey80",
                  force=30,
                  point.padding=1) +
  

  
  #create 2 geom_point plots, mapped from same data. 
  #Layering the 2 identical allows for customizing "outline" and "fill" of a given point. 
  #shape 16 is a solid circle. 
  #shape 1 is the outline of a circle.
  geom_point(shape=16, alpha=.5)+
  geom_point(shape=1)+
  

  
  scale_size(range = c(2, 15))+
  
  ggtitle("bubble ex")+

  theme_classic()
