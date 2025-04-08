library(ggplot2)
library(grid)

# Step 1: Get label buffer from a reference plot
get_label_buffer <- function(plot, units = "cm") {
  g <- ggplotGrob(plot)
  grid.newpage()
  grid.draw(g)
  
  # Get panel layout
  panel_ids <- grep("^panel", g$layout$name)
  panel_layout <- g$layout[panel_ids, ]
  
  panel_height_cm <- convertHeight(sum(g$heights[unique(panel_layout$t)]), units, valueOnly = TRUE)
  panel_width_cm  <- convertWidth(sum(g$widths[unique(panel_layout$l)]), units, valueOnly = TRUE)
  
  total_height_cm <- convertHeight(sum(g$heights), units, valueOnly = TRUE)
  total_width_cm  <- convertWidth(sum(g$widths), units, valueOnly = TRUE)
  
  list(
    label_buffer_height = total_height_cm - panel_height_cm,
    label_buffer_width  = total_width_cm - panel_width_cm
  )
}

# Step 2: Use the fixed label buffer in all aspect ratio calculations
gg_output_aspect_fixedbuffer <- function(plot,
                                         panel_width = 15,
                                         row_height = 0.4,
                                         n_rows,
                                         label_buffer,
                                         units = "cm",
                                         debug = FALSE) {
  total_height <- n_rows * row_height + label_buffer$label_buffer_height
  total_width  <- panel_width + label_buffer$label_buffer_width
  
  if (debug) {
    cat("Panel height:", n_rows * row_height,
        "| Label buffer height:", label_buffer$label_buffer_height, "\n")
    cat("Panel width :", panel_width,
        "| Label buffer width :", label_buffer$label_buffer_width, "\n")
    cat("Total size:", round(total_width, 2), "x", round(total_height, 2), units, "\n")
  }
  
  total_width / total_height
}





pad_labels <- function(label_vec, width = NULL) {
  if (is.null(width)) width <- max(nchar(label_vec))
  sprintf(paste0("%", width, "s"), label_vec)
}