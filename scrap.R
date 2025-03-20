





de2<-7

gost_v1[[de2]]$plot1 <- ggplot(data3, aes(x = jittered_x, 
                                          y = -log10(p_value), 
                                          label = g_label)) +
  
  # Add jittered points
  
  
  # Add jittered points
  geom_point(aes(size = intersection_size, 
                 color = pathway_list),
             shape = 16,
             alpha = 0.4)+# Ensures both fill & stroke work) +      
  geom_point(aes(size = intersection_size, 
                 color = pathway_list),
             shape = 1,
             stroke=.8,
             alpha = 0.8)+# Ensures both fill & stroke work) +  
  
  
  geom_label_repel(data=subset(data3),
                   aes(label = g_label,
                       point.size = intersection_size),
                   point.padding = .1,
                   #family="Open Sans",
                   alpha=.8,
                   size=3,
                   min.segment.length = 0,
                   segment.alpha=.25,
                   max.overlaps = 12,
                   force = 5,
                   hjust= 0.5,
                   vjust=0.5,
                   box.padding = 0,
                   show.legend = FALSE)+
  
  geom_hline(yintercept = 5, alpha=.6, linetype='dotted')+
  
  
  scale_x_continuous(
    breaks = breaks_dynamic,  # Breaks at each integer position
    limits = limits_dynamic,    # Padding on both ends
    labels = pathway_list_factor$pathway_list # Custom labels
  )+
  
  labs(x="pathway term group", 
       y="-log10(p value)",
       title = str_wrap(names(gost_v1)[[de2]], width = 30),
       intersection_size="Intersection Size")+
  
  guides(fill=guide_legend(title="Pathway Term Group"),
         size=guide_legend(title="Intersection Size"))+
  
  continuous_scale(
    aesthetics = c("size", "point.size"),  # Scale both point size and label size
    palette = my_pal(c(4, 12))  # Custom palette for scaling
  ) +
  coord_cartesian(clip = "off")+
  theme_classic()
  #theme(text = element_text(family = "Open Sans", size = 2, lineheight=8))


grid.newpage()
grid.draw(gost_v1[[7]]$plot1)
























# 
# 
# 
# 
# 
# 
# 
# gost_v1[[7]]$plot2<-grid.arrange(
#   gost_v1[[7]]$plot1, 
#   text_grob, 
#   ncol = 1, 
#   heights = c(3, 0.3),  # Maintain spacing
#   layout_matrix = rbind(c(1), c(2))
# )# Keeps layout structure stable

