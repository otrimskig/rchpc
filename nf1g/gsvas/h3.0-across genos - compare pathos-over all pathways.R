





# Ensure factor levels don’t interfere
df2 <- df2 %>%
  mutate(Pathway = as.character(Pathway))

# Step 1: Order by pathway_grouping and descending disease_mean
ordered_pathways <- df2 %>%
  distinct(Pathway, pathway_grouping, disease_mean) %>%
  arrange(pathway_grouping, desc(disease_mean)) %>%
  pull(Pathway)

# Step 2: Create chunks of 1000 based on the above order
pathway_chunks <- split(ordered_pathways, ceiling(seq_along(ordered_pathways) / 1000))

# Step 3: Output directory
dir.create("nf1g/gsvas/plots/chunks", recursive = TRUE, showWarnings = FALSE)

# Loop through and plot each chunk
for (i in seq_along(pathway_chunks)) {
  
  df_chunk <- df2 %>%
    filter(Pathway %in% pathway_chunks[[i]]) %>%
    mutate(Pathway = factor(Pathway, levels = rev(pathway_chunks[[i]])))  # reversed order for top-down
  
  p <- ggplot(data = df_chunk, aes(x = z_score, y = Pathway, color = patho_cat_name)) + 
    geom_point(aes(x = mean_z_score, y = Pathway), color = "black", shape = 124, size = 3) +
    geom_point(size = 2, alpha = 0.5) +
    ggh4x::facet_nested(pathway_grouping ~ patho_cat_name, scales = "free_y", space = "free") +
    theme_bw() +
    theme(
      strip.placement = "outside",
      strip.text = element_text(size = 10),
      strip.background = element_blank(),
      axis.text.y = element_text(size = 7.5),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      panel.spacing = unit(0.0, "lines"),
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    scale_y_discrete(expand = expansion(add = c(1, 1))) +
    scale_color_manual(values = col_map[["hist_cat_name"]])
  
  # Save plot
  ggsave(
    filename = sprintf("nf1g/gsvas/plots/chunks/plot_chunk_%02d.pdf", i),
    plot = p,
    width = 2,
    height = 14,
    scale = 8,
    limitsize = FALSE
  )
  
  cat(cli::col_blue(i), "of", cli::col_red(length(pathway_chunks)), "\n")
}

















plot_paths <- fs::dir_info("nf1g/gsvas/plots/chunks", regexp = "\\.pdf$", recurse = FALSE) %>%
  filter(modification_time>Sys.time()-lubridate::minutes(10))%>%
  filter(type == "file") %>%
  arrange(path) %>%
  tibble()



qpdf::pdf_combine(plot_paths$path, output = "nf1g/gsvas/plots/chunks/combined_chunks.pdf")
