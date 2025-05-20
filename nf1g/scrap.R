# # Step 1: Prep the data
# library(ggrepel)
# library(ggplot2)
# 
# label_data <- subset(data3, !is.na(g_label))
# 
# # Step 2: Compute repelled positions manually
# repel_positions <- ggrepel:::calc_repel_positions(
#   data = label_data,
#   force = 5,
#   point.padding = 0.1,
#   box.padding = 0,
#   segment.length = 0,
#   max.iter = 2000,
#   direction = "both",
#   xlim = range(data3$your_x_var),
#   ylim = range(data3$your_y_var),
#   x.name = "your_x_var",
#   y.name = "your_y_var"
# )
# 
# # Step 3: Add the computed x/y to your data
# label_data$x_rep <- repel_positions$x
# label_data$y_rep <- repel_positions$y
# 
# # Step 4: Plot — now use geom_label() twice, perfectly aligned
# ggplot(data3, aes(x = your_x_var, y = your_y_var)) +
#   geom_point() +  # your base plot layer, optional
#   
#   # Background label fill with low alpha
#   geom_label(data = label_data,
#              aes(x = x_rep, y = y_rep, label = g_label, fill = pathway_list),
#              color = NA,
#              alpha = 0.2,
#              size = 3,
#              show.legend = FALSE) +
#   
#   # Foreground text + outline with higher alpha
#   geom_label(data = label_data,
#              aes(x = x_rep, y = y_rep, label = g_label, fill = pathway_list),
#              color = "black",
#              alpha = 0.7,
#              size = 3,
#              show.legend = FALSE)
# 










repel_boxes <- getFromNamespace("repel_boxes", "ggrepel")


calc_repel_positions <- function(data,
                                 force = 1,
                                 point.padding = 1e-6,
                                 box.padding = 0.25,
                                 segment.length = 0.1,
                                 max.iter = 2000,
                                 direction = "both",
                                 xlim = c(-Inf, Inf),
                                 ylim = c(-Inf, Inf),
                                 x.name = "x",
                                 y.name = "y") {
  
  if (!requireNamespace("grid")) {
    stop("This function requires the 'grid' package.")
  }
  
  if (!requireNamespace("ggrepel")) {
    stop("This function depends on 'ggrepel'.")
  }
  
  # minimal pre-processing
  data$x <- data[[x.name]]
  data$y <- data[[y.name]]
  data$label <- as.character(data$label)
  
  repel_df <- ggrepel:::repel_boxes(
    data$x, data$y, data$label,
    force = force,
    point.padding = point.padding,
    box.padding = box.padding,
    segment.length = segment.length,
    max.iter = max.iter,
    direction = direction,
    xlim = xlim,
    ylim = ylim
  )
  
  return(data.frame(x = repel_df$x1, y = repel_df$y1))
}
