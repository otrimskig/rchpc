library(grid)
library(gridExtra)
library(ggplot2)

gg_output_aspect <- function(plot,
                             panel_width = 15,
                             row_height = 0.4,
                             n_rows,
                             units = "cm",
                             debug = FALSE) {
  # Build grob and draw it to resolve layout
  g <- ggplotGrob(plot)
  grid.newpage()
  grid.draw(g)
  
  # Get panel layout info
  panel_ids <- grep("^panel", g$layout$name)
  panel_layout <- g$layout[panel_ids, ]
  
  # Sum up actual panel heights
  panel_heights <- g$heights[unique(panel_layout$t)]
  panel_height_cm <- convertHeight(sum(panel_heights), units, valueOnly = TRUE)
  total_height_cm <- convertHeight(sum(g$heights), units, valueOnly = TRUE)
  label_buffer_height <- total_height_cm - panel_height_cm
  
  # Same for width (in case of nested strips or y labels)
  panel_widths <- g$widths[unique(panel_layout$l)]
  panel_width_cm <- convertWidth(sum(panel_widths), units, valueOnly = TRUE)
  total_width_cm <- convertWidth(sum(g$widths), units, valueOnly = TRUE)
  label_buffer_width <- total_width_cm - panel_width_cm
  
  # Use user-specified desired panel size and row count
  desired_panel_height <- n_rows * row_height
  total_height <- desired_panel_height + label_buffer_height
  total_width <- panel_width + label_buffer_width
  
  if (debug) {
    cat("Panel height:", desired_panel_height, "Label buffer height:", label_buffer_height, "\n")
    cat("Panel width :", panel_width, "Label buffer width :", label_buffer_width, "\n")
    cat("Total size:", round(total_width, 2), "x", round(total_height, 2), units, "\n")
  }
  
  # Return the aspect ratio to use in saving
  total_width / total_height
}
