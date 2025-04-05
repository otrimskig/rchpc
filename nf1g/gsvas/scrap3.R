

# 
# df2<-df1
# 
# ggplot(data=df1, aes(x=z_score, y=reorder(Pathway, disease_mean)))+
#   geom_point(aes(color=patho_cat_name), size=3, alpha=.5, position = position_dodge(width = 1, preserve = "single"))+
#   theme_bw()
#   #facet_grid(cols = vars(patho_cat_name)) 
# ggplot(data = df1, aes(x = z_score, y = reorder(Pathway, disease_mean))) +
#   geom_point(aes(color = patho_cat_name), size = 3, alpha = 0.5) +
#   theme_bw() +
#   facet_wrap(~ patho_cat_name, scales = "free_y", ncol = 1)
# 


# 
# ggplot(data = df1, aes(x = z_score, y = interaction(Pathway, patho_cat_name), color = patho_cat_name)) +
#   geom_point(size = 3, alpha = 0.5) +
#   theme_bw()
# ggplot(df1, aes(x = z_score, y = Pathway, color = patho_cat_name)) +
#   geom_point(size = 3, alpha = 0.5, position = "identity") +  # No dodging, stacked points
#   theme_bw() +
#   theme(axis.text.y = element_text(size = 8))  # Optional: Adjust y-axis text size for clarity
# 
# 
# 






# 
# ggplot(data = df1, aes(x = z_score, y = patho_cat_name, color = patho_cat_name)) + 
#   geom_point(position = position_dodge(width = 1), size = 3, alpha = 0.6) +
#   theme_bw() +
#   theme(
#     axis.text.y = element_text(angle = 0, hjust = 1),  # Adjust y-axis text alignment
#     panel.grid.major.y = element_line(color = "gray", size = 0.2),  # Major gridlines
#     panel.grid.minor.y = element_line(color = "lightgray", size = 0.5),  # Minor gridlines
#     panel.grid.major.x = element_blank(),  # Remove major gridlines on x-axis
# )+ facet_wrap(~ Pathway, ncol = 1)
#   
# 
# 
# 
# 







# 
# 
# 
# ggplot(data = df1, aes(x = z_score, y = Pathway, color = patho_cat_name)) + 
#   geom_point(position = position_dodge(width = 1), size = 3, alpha = 0.6) +
#   theme_bw() +
#   theme(
#     axis.text.y = element_text(angle = 0, hjust = 1),  # Adjust y-axis text alignment
#     panel.grid.major.x = element_blank(),  # Remove gridlines on x-axis
#     panel.grid.minor.x = element_blank(),  # Remove minor gridlines on x-axis
#     panel.grid.major.y = element_blank(),  # Remove default gridlines on y-axis
#     panel.grid.minor.y = element_blank()   # Remove minor gridlines
#   ) +
#   scale_y_discrete(expand = c(0, 0)) +  # Ensure no extra space at top and bottom
#   # Add horizontal lines between pathways
#   geom_tile(aes(x = -Inf, xend = Inf, y = Pathway, yend = Pathway), 
#             color = "gray", size = 0.5, inherit.aes = FALSE)
# 
# 
# ggplot(data = df1, aes(x = z_score, 
#                        y = reorder(Pathway, disease_mean), 
#                        color = patho_cat_name)) +
#   geom_point(size = 3, alpha = 0.5, 
#              position = position_dodge2(width = 0.1, preserve = "single")) +
#   theme_bw()
# 
# 




# 
# 
# 
# ggplot(data = df1, aes(x = z_score, y = Pathway, color = patho_cat_name)) + 
#   geom_point(position = position_dodge(width = .75), size = 2, alpha = 0.6) +
#   theme_bw() +
#   theme(
#     axis.text.y = element_text(angle = 0, hjust = 1),  # Adjust y-axis text alignment
#     panel.grid.major.x = element_blank(),  # Remove gridlines on x-axis
#     panel.grid.minor.x = element_blank(),  # Remove minor gridlines on x-axis
#     panel.grid.major.y = element_blank(),  # Remove default gridlines on y-axis
#     panel.grid.minor.y = element_blank()   # Remove minor gridlines
#   ) +
#   scale_y_discrete(minor.breaks = 1.5:6.5) 
#   #scale_y_discrete(expand = c(0, 0)) +  # Ensure no extra space at top and bottom
#   # Add horizontal lines between pathways


# 
# # Create a numeric mapping of y-axis categories for horizontal lines
# y_levels <- levels(factor(df1$patho_cat_name))
# y_positions <- seq_along(y_levels)
# 
# # Make a data frame for horizontal lines
# hline_data <- data.frame(
#   y = y_positions
# )
# 
# ggplot(data = df1, aes(x = z_score, y = patho_cat_name, color = patho_cat_name)) + 
#   # Horizontal lines at each y-category
#   geom_hline(data = hline_data, aes(yintercept = y), color = "grey85", linetype = "solid", alpha=.5) +
#   
#   # The main points
#   geom_point(position = position_dodge(width = 1), size = 1, alpha = 0.5) +
#   
#   facet_grid(Pathway ~ ., switch = "y") +
#   theme_bw() +
#   theme(
#     strip.placement = "outside",
#     strip.text.y.left = element_text(angle = 0, hjust = 1),
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     panel.spacing = unit(0, "lines"),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_line(color = "grey90", linetype = "dotted"),
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank()
#   ) +
#   scale_y_discrete(expand = expansion(add = 1.5))  # Some breathing room





# ggplot(data = df1, aes(x = z_score, y = 0, color = patho_cat_name)) + 
#   geom_point(size = 2, alpha = 0.6) +
#   facet_grid(Pathway ~ patho_cat_name, switch = "y") +
#   theme_bw() +
#   theme(
#     strip.placement = "outside",
#     strip.text.y.left = element_text(angle = 0, hjust = 1),
#     axis.text.y = element_blank(),
#     axis.ticks.y = element_blank(),
#     axis.title.y = element_blank(),
#     panel.spacing = unit(0, "lines"),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     plot.margin = margin(5, 5, 5, 5),                   # Optional: tighter outer margins
#     strip.text = element_text(size = 8),                # Smaller facet label text
#     strip.background = element_blank()                  # Optional: no background
#   ) +
#   scale_y_continuous(limits = c(-0.2, 0.2), expand = expansion(mult = 0.1)) +
#   guides(color = "none")
# 
# 
