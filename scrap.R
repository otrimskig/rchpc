
# Ensure sample names match matrix column names
valid_a_samples <- intersect(group_a_samples, colnames(sub_matrix))
valid_b_samples <- intersect(group_b_samples, colnames(sub_matrix))

# Extract groups as matrices
group_a_matrix <- sub_matrix[, valid_a_samples, drop = FALSE]
group_b_matrix <- sub_matrix[, valid_b_samples, drop = FALSE]

# Compute mean per group
mean_a <- rowMeans(group_a_matrix, na.rm = TRUE)
mean_b <- rowMeans(group_b_matrix, na.rm = TRUE)

# Compute fold change
fold_changes <- mean_a - mean_b

# Compute standard deviation per group
std_a <- apply(group_a_matrix, 1, sd, na.rm = TRUE)
std_b <- apply(group_b_matrix, 1, sd, na.rm = TRUE)

# Apply t-test across all pathways (rows)
p_values <- apply(sub_matrix, 1, function(row) {
  a_values <- row[valid_a_samples]
  b_values <- row[valid_b_samples]
  
  # Ensure at least 2 non-NA values in both groups
  if (sum(!is.na(a_values)) > 1 & sum(!is.na(b_values)) > 1) {
    t.test(a_values, b_values)$p.value
  } else {
    NA  # Return NA if not enough data
  }
})

# Create results tibble
results <- tibble(
  Pathway = rownames(sub_matrix),
  p_value = p_values,
  mean_group_a = mean_a,
  mean_group_b = mean_b,
  std_group_a = std_a,
  std_group_b = std_b
) %>%
  arrange(p_value)  # Sort by significance

# Create the output list with metadata
analysis_output <- list(
  results = results,
  metadata = list(
    date = Sys.time(),
    group_comparison = paste(group_factors[1], "vs", group_factors[2]),
    group_a = paste(group_factors[1]),
    group_b = paste(group_factors[2]),
    group_a_samples = valid_a_samples,
    group_b_samples = valid_b_samples
  )
)




