source("libs.R")

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(gprofiler2)
library(ggpp)

data1<-readRDS("nf1g/gresult1test.rds")


data2 <- data1%>%
  mutate(point_index=1:n())%>%
  rename(pathway_list=source)



centered_jitter <- function(x, width = 0.4, seed = NULL) {
  if (!is.null(seed)) set.seed(seed) # Ensure reproducibility
  n <- length(x)
  if (n == 1) {
    return(x)  # Keep single points unchanged
  } else {
    jitter_offsets <- runif(n, min = -width / 2, max = width / 2)
    jitter_offsets <- jitter_offsets - mean(jitter_offsets)  # Center around 0
    return(x + jitter_offsets)
  }
}




# Apply to dataset
data3 <- data2 %>%
  mutate(group_x_coord= as.numeric(factor(pathway_list)))%>%
  group_by(pathway_list)%>%
  mutate(jittered_x = centered_jitter(group_x_coord, width = 0.4, seed = 3))%>%
  mutate(dist_from_group_x=jittered_x-group_x_coord)%>%
  
  mutate(label_x_pos=if_else(dist_from_group_x<0, group_x_coord-.3, group_x_coord+.3))%>%
  
  mutate(g_label = if_else(-log10(p_value) > 5|grepl("MHC", term_name), 
                           stringr::str_wrap(term_name, width=25),
                           ""))





my_pal <- function(range = c(2, 12)) {
  force(range)
  function(x) scales::rescale(x, to = range, from = c(0, 1))
}

x_labels<-data3%>%
  arrange(group_x_coord)%>%
  select(pathway_list)%>%distinct()%>%pull()





p2 <- ggplot(data3, aes(x = jittered_x, 
                        y = -log10(p_value), 
                        label = g_label)) +
  
  
  
  
  # Add jittered points
  geom_point(aes(size = intersection_size, 
                 fill = pathway_list),
              color = "black",
              alpha = 0.6, 
              shape = 21,    # Ensures both fill & stroke work
              stroke = 0.5) +  # Controls outline thickness
  
  #facet_wrap(vars(pathway_list), nrow = 1)+

  # geom_text_s(aes(label = g_label),
  #             position = position_nudge_to(x = c(1.1, 0.9)),
  #             box.padding = 0.5) +
  # 
  geom_label_repel(data=subset(data3),
                   aes(label = g_label,
                       point.size = intersection_size),
                       #color=pathway_list),
                   xlim=c(0,6),
                   #direction = "both",
                   point.padding = 0,
                   min.segment.length = 0,
                   #size = 3,
                   # segment.ncp = 3,
                   # segment.curvature = -.1,
                   max.overlaps = 9,
                   force = 7,
                   hjust= 0,
                   alpha=.3,
                   box.padding = 0,
                   show.legend = FALSE)+
  

  scale_x_continuous(breaks = seq(1, 5, by = 1),
                     limits = c(0,5.5),
                     labels=x_labels)+
  
  
  
  labs(x="pathway term group", 
      y="-log10(p value)",
      title = "4KO vs. ATRX wt GO Analysis",
      intersection_size="Intersection Size")+

  guides(fill=guide_legend(title="Pathway Term Group"),
         size=guide_legend(title="Intersection Size"))+
  
  continuous_scale(
    aesthetics = c("size", "point.size"),  # Scale both point size and label size
    palette = my_pal(c(2, 12))  # Custom palette for scaling
  ) +
  coord_cartesian(clip = "off")+
  theme_classic()
  


p2

# print(p2)


#save plot
plot_object_name_for_session<-"p2"
plot_filename_output<-"nf1g/gost/plots/atrx_comp_gost.pdf"


metadata_text <- paste0("src: ", 
                        rstudioapi::getSourceEditorContext()$path %>%
                          sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.) %>%
                          sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.), 
                        "\n", "time: ", 
                        lubridate::round_date(Sys.time(), "second"))

text_grob <- grid::textGrob(metadata_text, x=.05, just="left", gp = grid::gpar(fontsize = 9, col = "gray30"))


p3src<-gridExtra::grid.arrange(get(plot_object_name_for_session), text_grob, ncol = 1, heights = c(3, 0.3))

ggsave(plot_filename_output,
       
       plot=p3src,
       limitsize = FALSE,
       
       height=5,
       width=8,
       scale = 1.1,
       dpi=600)

