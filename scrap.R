

manual_legend<-

spec_df1 %>%
  #arrange(p_value) %>%
  mutate(group_label=resultant_geno)%>%
  group_by(group_label)%>%
  count(group_label)%>%
  ungroup()%>%
  
  
  mutate(
    group_a = factor(group_label, levels = unique(group_label[order(group_label)])),
    #group_b = factor(group_b, levels = unique(group_b[order(p_value)])),
    row_id = row_number()
  )%>%
  
  ggplot() +
  
  #set up tiles for comparisons.
  #set up top/bottom color comparison.
  geom_tile(aes(x = 0, y = 0), fill = "#ccffff", alpha = .0, width = 10, height = nrows) +  # Empty tile
  geom_tile(aes(x = -5, y = nrows-row_id, fill = group_a), width = 1.25, height = 0.25, color = "#5c5c5c") +
  #geom_tile(aes(x = -5, y = nrows-row_id-.2, fill = group_b), width = 1.25, height = 0.25, color = "#5c5c5c") +
  
  
  #text for each comparison.
  geom_text(aes(
    x = -4, 
    y = nrows-row_id, 
    
    label = group_a, 
    ), 
    color =  "black", 
    size = 7, hjust = 0) +

  scale_color_identity()+
  
  scale_fill_manual(values = col_map[[split_curves_by]]) +
  theme_minimal() +
  
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0))+
  coord_fixed(ratio = 1)
