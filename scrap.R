df <- data3 %>%
  group_by(pathway_list) %>%
  mutate(
    rank = rank(p_value, ties.method = "average"),  # Rank p-values within each pathway
    percentile = (1 - percent_rank(p_value)) * 100  # Convert to percentile (0-100 scale)
  ) %>%
  ungroup()



plot_obj <- gost_v1[["dexp-patho_cat_name-Spindle and epithelioid tumors-nf1 KO; pten KO; ink KO; atrx KO v. nf1 KO; pten KO; ink KO; atrx wt.rds"]][["plot1"]]

# Extract colors from the ggplot object
plot_colors <- ggplot_build(plot_obj)$data[[1]]$colour  

# Show unique colors
unique(plot_colors)
