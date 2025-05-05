library(ggplot2)
library(patchwork)  # for plot layout

# Ensure the order of genotypes (if needed)
df_props_stack1$resultant_geno <- factor(df_props_stack1$resultant_geno)
df_stack_pen$resultant_geno <- factor(df_stack_pen$resultant_geno, 
                                      levels = levels(df_props_stack1$resultant_geno))

pstack2 <- ggplot(df_props_stack1, aes(y = resultant_geno, 
                                       x = perc, 
                                       fill = hist_grade_name)) +
  
  # Main stacked bars, now horizontal by default
  geom_bar(stat = "identity", position = "stack", width=.7, alpha = 0.7) +
  
  scale_fill_manual(values = col_map$hist_grade_name, name = NULL)+
  
  # Tile annotations now placed at x = pen
  geom_tile(
    data = df_stack_pen,
    aes(y = resultant_geno, x = pen),
    inherit.aes = FALSE,
    fill = "black",
    height = 0.9,
    width = 0.2,
    alpha = 1
  ) +
  

  # White label background
  geom_label(
    data = df_stack_pen,
    aes(y = resultant_geno, x = pen+2.2, label = paste0(round(pen), "%")),
    inherit.aes = FALSE,
    fill = "white",
    color = NA,
    label.r = unit(0, "pt"),
    size = 3
  ) +
  scale_y_discrete(limits = geno_order)+
  
  
  # Black overlay text
  geom_text(
    data = df_stack_pen,
    aes(y = resultant_geno, x = pen+2.2, label = paste0(round(pen), "%")),
    inherit.aes = FALSE,
    color = "black",
    size = 3
  ) +
  
  
  
  
  
  
  
  theme_pubr() +
  theme(
    axis.text.y =  element_blank(),
    plot.margin = margin(5, 50, 5, 0),
    axis.text.x = element_text(size = 12),
    guides(fill = guide_legend(nrow = 9))
  ) +
  
  labs(x = "% of cohort", y = NULL) +
  
  scale_x_continuous(expand = expansion(mult = c(0, 0)))
  
  
pstack2
  

geno_order <- levels(df_props_stack1$resultant_geno)

df_anno <- data.frame(
  resultant_geno = factor(geno_order, levels = geno_order),  # Explicitly factor with correct order
  y = geno_order,
  x = 1,
  color = df_stack_pen$hex_col  # assuming it's in the same order
)


# Annotation plot (one row per genotype)
p_anno <- ggplot(df_anno, aes(y = resultant_geno)) +

  geom_tile(aes(x = 9.9, fill = color), width = .2, height = 0.8) +

  geom_text(aes(x = 9.6, label = resultant_geno), hjust = 1, size = 4) +
  
  scale_fill_identity() +
  scale_y_discrete(limits = geno_order)+
  #coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    #plot.margin = margin(5, 10, 5, 5),
    axis.text = element_blank()
  ) +
 xlim(0, 10)


combined_plot <- p_anno + pstack2 + plot_layout(widths = c(1, 3))


combined_plot

