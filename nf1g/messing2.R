


comps %>%
  arrange(p_value) %>%
  mutate(
    group_a = factor(group_a, levels = unique(group_a[order(p_value)])),
    group_b = factor(group_b, levels = unique(group_b[order(p_value)])),
    row_id = row_number()
  )%>%
  
  ggplot() +
  
  geom_tile(aes(x = 7.5, y = 7.5), fill = "#ccffff", alpha = .50, width = 15, height = 15) +  # Empty tile
  
  
  #set sizes of legend keys.
  geom_tile(aes(x = 0.75, y = row_id+.15, fill = group_a), 
            width = 1.5, height = 0.15, color = "#5c5c5c") +
  geom_tile(aes(x = .75, y = row_id-.15, fill = group_b), 
            width = 1.5, height = 0.15, color = "#5c5c5c") +
  
  
  geom_text(aes(
    x = 2, y = row_id, 
    
    label = ifelse(p_value < 0.001, paste0("***P", 
                                           format(p_value, scientific = FALSE, digits = 4))
      
  
                   
                   format(p_value, scientific = F, digits = 4)),
      
      
      
    fontface = ifelse(p_value <= 0.05, "bold", "plain"),
    
    color = ifelse(p_value <= 0.05, "black", "gray50")
  ), size = 10, hjust = 0) +
  
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
