source("libs.R")
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
  
  # Add plot margin handling
  plot_margin <- g$theme$plot.margin
  
  # Ensure margins are numeric values, if they're not unit objects
  if (is.null(plot_margin)) {
    plot_margin <- unit(c(0, 0, 0, 0), "cm")  # Default margin if none defined
  }
  
  # Extract the margin values and convert them to the requested unit
  margin_top <- convertUnit(plot_margin[1], units, valueOnly = TRUE)
  margin_right <- convertUnit(plot_margin[2], units, valueOnly = TRUE)
  margin_bottom <- convertUnit(plot_margin[3], units, valueOnly = TRUE)
  margin_left <- convertUnit(plot_margin[4], units, valueOnly = TRUE)
  
  # Adjust total width and height to account for margins
  total_width <- panel_width + label_buffer_width + margin_left + margin_right
  total_height <- n_rows * row_height + label_buffer_height + margin_top + margin_bottom
  
  if (debug) {
    cat("Panel height:", n_rows * row_height, "Label buffer height:", label_buffer_height, "\n")
    cat("Panel width :", panel_width, "Label buffer width :", label_buffer_width, "\n")
    cat("Margins:", "Top =", margin_top, "Right =", margin_right, "Bottom =", margin_bottom, "Left =", margin_left, "\n")
    cat("Total size:", round(total_width, 2), "x", round(total_height, 2), units, "\n")
  }
  
  # Return the aspect ratio to use in saving
  total_width / total_height
}






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



