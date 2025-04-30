library(nnet)
df_props1$hist_grade_name <- relevel(df_props1$hist_grade_name, ref = "Grade 1")

model <- multinom(hist_grade_name ~ resultant_geno, weights = n, data = df_props1)
summary(model)


# Predict probabilities
new_data <- data.frame(resultant_geno = levels(df_props1$resultant_geno))
predicted_probs <- predict(model, newdata = new_data, type = "probs")

# Convert to long format for ggplot
library(tidyr)
library(dplyr)

pred_df <- as.data.frame(predicted_probs)
pred_df$resultant_geno <- new_data$resultant_geno

pred_long <- pivot_longer(pred_df, 
                          cols = -resultant_geno, 
                          names_to = "hist_grade_name", 
                          values_to = "probability")

# Plot
library(ggplot2)

ggplot(pred_long, aes(x = resultant_geno, y = probability, fill = hist_grade_name)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Predicted Distribution of Histological Grades by Genotype",
       x = "Genotype", y = "Predicted Proportion", fill = "Hist Grade") +
  theme_minimal() +
  coord_flip()


# Odds ratios
exp(coef(model))


coefs <- summary(model)$coefficients
ses <- summary(model)$standard.errors
z_scores <- coefs / ses
p_values <- 2 * (1 - pnorm(abs(z_scores)))
p_values






library(nnet)
library(dplyr)

# Get all unique genotype levels
geno_levels <- levels(df_props1$resultant_geno)

# Create all unique pairs
geno_pairs <- combn(geno_levels, 2, simplify = FALSE)

# Function to run multinom and return p-value
pairwise_test <- function(pair, data) {
  sub_data <- filter(data, resultant_geno %in% pair)
  sub_data <- droplevels(sub_data)
  
  full_model <- multinom(hist_grade_name ~ resultant_geno, weights = n, data = sub_data, trace = FALSE)
  null_model <- multinom(hist_grade_name ~ 1, weights = n, data = sub_data, trace = FALSE)
  
  LRT <- 2 * (logLik(full_model) - logLik(null_model))
  df <- attr(logLik(full_model), "df") - attr(logLik(null_model), "df")
  p_val <- pchisq(LRT, df = df, lower.tail = FALSE)
  
  data.frame(geno1 = pair[1], geno2 = pair[2], p_value = p_val)
}

# Run for all pairs
results <- do.call(rbind, lapply(geno_pairs, pairwise_test, data = df_props1))

# Adjust for multiple comparisons
results$p_adj <- p.adjust(results$p_value, method = "BH")

results %>% arrange(geno1, geno2)


















library(dplyr)
library(tidyr)
library(purrr)

# Get unique genotypes and grades
genos <- levels(df_props1$resultant_geno)
grades <- levels(df_props1$hist_grade_name)

# Function to run Fisher's test for a single grade/genotype pair
fisher_for_pair <- function(grade, geno_pair, data) {
  g1 <- geno_pair[1]
  g2 <- geno_pair[2]
  
  sub_data <- data %>% 
    filter(resultant_geno %in% c(g1, g2)) %>%
    group_by(resultant_geno, hist_grade_name) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    complete(resultant_geno = c(g1, g2), hist_grade_name = grades, fill = list(n = 0)) %>%
    mutate(grade_of_interest = ifelse(hist_grade_name == grade, "target", "other")) %>%
    group_by(resultant_geno, grade_of_interest) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    pivot_wider(names_from = grade_of_interest, values_from = n, values_fill = 0) %>%
    column_to_rownames("resultant_geno") %>%
    as.matrix()
  
  if (all(rowSums(sub_data) > 0) && all(colSums(sub_data) > 0)) {
    test <- fisher.test(sub_data)
    data.frame(grade = grade, geno1 = g1, geno2 = g2, 
               p_value = test$p.value, odds_ratio = test$estimate)
  } else {
    data.frame(grade = grade, geno1 = g1, geno2 = g2, 
               p_value = NA, odds_ratio = NA)
  }
}

# Generate all genotype pairs
geno_pairs <- combn(genos, 2, simplify = FALSE)

# Run the tests
results_fisher <- map_dfr(grades, function(g) {
  map_dfr(geno_pairs, function(p) fisher_for_pair(g, p, df_props1))
})

# Adjust p-values within each grade for multiple comparisons
results_fisher <- results_fisher %>%
  group_by(grade) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(grade, p_adj)

results_fisher




library(ggplot2)
library(forcats)

# Create a label for genotype pair
results_fisher <- results_fisher %>%
  mutate(pair_label = paste(geno1, "vs", geno2),
         log_odds = log2(odds_ratio),
         signif = case_when(
           is.na(p_adj) ~ "NA",
           p_adj < 0.01 ~ "***",
           p_adj < 0.05 ~ "**",
           p_adj < 0.1 ~ "*",
           TRUE ~ ""
         ))

# Plot heatmap of log2(odds ratio), annotated with significance
ggplot(results_fisher, aes(x = grade, y = fct_rev(pair_label), fill = log_odds)) +
  geom_tile(color = "white") +
  geom_text(aes(label = signif), size = 3) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                       midpoint = 0, na.value = "grey80",
                       name = "log2(OR)") +
  labs(title = "Pairwise Fisher's Tests by Grade and Genotype",
       x = "Histological Grade", y = "Genotype Comparison") +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1))












library(ggplot2)
library(forcats)

# Prepare data
plot_data <- results_fisher %>%
  filter(!is.na(odds_ratio)) %>%
  mutate(
    pair = paste(geno1, "vs", geno2),
    pair = fct_reorder(pair, log2(odds_ratio)),
    log2_or = log2(odds_ratio),
    sig_level = cut(p_adj,
                    breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
                    labels = c("***", "**", "*", "ns"))
  )

# Plot
ggplot(plot_data, aes(x = log2_or, y = pair, color = sig_level)) +
  geom_point(aes(size = -log10(p_adj))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~grade, scales = "free_y") +
  scale_color_manual(values = c("***" = "#d73027", "**" = "#fc8d59", "*" = "#fee08b", "ns" = "grey70")) +
  scale_size_continuous(range = c(1, 6)) +
  labs(
    title = "Pairwise Comparison of Histological Grades Between Genotypes",
    x = "log2(Odds Ratio) (positive = more in geno1)",
    y = "Genotype Pair",
    color = "Significance",
    size = "-log10(adj p)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")










library(ggplot2)
library(forcats)

# Clean and prep data
plot_data_A <- results_fisher %>%
  filter(!is.na(odds_ratio)) %>%
  mutate(
    pair = paste(geno1, "vs", geno2),
    log2_or = log2(odds_ratio),
    sig_level = cut(p_adj,
                    breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
                    labels = c("***", "**", "*", "ns"))
  )

ggplot(plot_data_A, aes(x = grade, y = log2_or, color = sig_level, size = -log10(p_adj))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_point() +
  facet_wrap(~pair, scales = "free_x") +
  scale_color_manual(values = c("***" = "#d73027", "**" = "#fc8d59", "*" = "#fee08b", "ns" = "gray80")) +
  scale_size_continuous(range = c(1.5, 6)) +
  labs(
    title = "Plot A: Grade Differences by Genotype Pair",
    x = "Histological Grade",
    y = "log2(Odds Ratio)",
    color = "Significance",
    size = "-log10(adj p)"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Ensure 'resultant_geno' is included in the results_fisher dataframe
plot_data_grid <- results_fisher %>%
  filter(!is.na(odds_ratio)) %>%
  mutate(
    log2_or = log2(odds_ratio),
    sig_level = cut(p_adj,
                    breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
                    labels = c("***", "**", "*", "ns")),
    geno1 = factor(geno1, levels = sort(unique(df_props1$resultant_geno))),
    geno2 = factor(geno2, levels = sort(unique(df_props1$resultant_geno)))
  ) %>%
  filter(as.integer(geno1) < as.integer(geno2))  # only upper triangle to avoid dupes


# Plot
ggplot(plot_data_grid, aes(x = grade, y = log2_or, color = sig_level, size = -log10(p_adj))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_point() +
  facet_grid(rows = vars(geno1), cols = vars(geno2)) +
  scale_color_manual(values = c("***" = "#d73027", "**" = "#fc8d59", "*" = "#fee08b", "ns" = "gray80")) +
  scale_size_continuous(range = c(1.5, 6)) +
  labs(
    title = "Plot A (Improved): Genotype-by-Genotype Matrix of Grade Differences",
    x = "Histological Grade",
    y = "log2(Odds Ratio)",
    color = "Significance",
    size = "-log10(adj p)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



















# Clean and prep data
plot_data_B <- plot_data_A %>%
  mutate(pair = fct_reorder(pair, log2_or))  # Optional: reorder within facets

ggplot(plot_data_B, aes(x = pair, y = log2_or, color = sig_level, size = -log10(p_adj))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_point() +
  facet_wrap(~grade, scales = "free_x") +
  scale_color_manual(values = c("***" = "#d73027", "**" = "#fc8d59", "*" = "#fee08b", "ns" = "gray80")) +
  scale_size_continuous(range = c(1.5, 6)) +
  labs(
    title = "Plot B: Genotype Pair Differences by Grade",
    x = "Genotype Pair",
    y = "log2(Odds Ratio)",
    color = "Significance",
    size = "-log10(adj p)"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))







# Prepare data for Plot B (facet by grade, showing all genotype comparisons)
plot_data_B <- results_fisher %>%
  filter(!is.na(odds_ratio)) %>%
  mutate(
    log2_or = log2(odds_ratio),
    sig_level = cut(p_adj,
                    breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
                    labels = c("***", "**", "*", "ns")),
    geno1 = factor(geno1, levels = sort(unique(df_props1$resultant_geno))),
    geno2 = factor(geno2, levels = sort(unique(df_props1$resultant_geno))),
    pair = paste(geno1, "vs", geno2)  # Creating the pair variable for comparisons
  ) %>%
  filter(as.integer(geno1) < as.integer(geno2))  # only upper triangle to avoid duplicate comparisons

# Plot for facet by grade
ggplot(plot_data_B, aes(x = pair, y = log2_or, color = sig_level, size = -log10(p_adj))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_point() +
  facet_wrap(~grade, scales = "free_x") +
  scale_color_manual(values = c("***" = "#d73027", "**" = "#fc8d59", "*" = "#fee08b", "ns" = "gray80")) +
  scale_size_continuous(range = c(1.5, 6)) +
  labs(
    title = "Plot B: Genotype Pair Differences by Grade",
    x = "Genotype Pair",
    y = "log2(Odds Ratio)",
    color = "Significance",
    size = "-log10(adj p)"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





# Prepare data for Plot B (facet by grade, showing all genotype comparisons in a matrix)
plot_data_B <- results_fisher %>%
  filter(!is.na(odds_ratio)) %>%
  mutate(
    log2_or = log2(odds_ratio),
    sig_level = cut(p_adj,
                    breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
                    labels = c("***", "**", "*", "ns")),
    geno1 = factor(geno1, levels = sort(unique(df_props1$resultant_geno))),
    geno2 = factor(geno2, levels = sort(unique(df_props1$resultant_geno))),
    pair = paste(geno1, "vs", geno2)  # Create pair variable for plotting
  ) %>%
  filter(as.integer(geno1) < as.integer(geno2))  # only upper triangle to avoid duplicate comparisons

# Plot for facet by grade with matrix-style layout (genotype pair matrix within each grade)
ggplot(plot_data_B, aes(x = geno2, y = geno1, fill = log2_or)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", log2_or)), color = "white", size = 4) +
  facet_wrap(~grade, scales = "free") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  labs(
    title = "Plot B: Genotype Pair Differences by Grade (Matrix Layout)",
    x = "Genotype (geno2)",
    y = "Genotype (geno1)",
    fill = "log2(Odds Ratio)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12)
  )
















# Prepare data for Plot B (facet by grade, showing all genotype comparisons in a matrix)
plot_data_B <- results_fisher %>%
  filter(!is.na(odds_ratio)) %>%
  mutate(
    log2_or = log2(odds_ratio),
    sig_level = cut(p_adj,
                    breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
                    labels = c("***", "**", "*", "ns")),
    geno1 = factor(geno1, levels = sort(unique(df_props1$resultant_geno))),
    geno2 = factor(geno2, levels = sort(unique(df_props1$resultant_geno))),
    pair = paste(geno1, "vs", geno2)  # Create pair variable for plotting
  ) %>%
  filter(as.integer(geno1) < as.integer(geno2))  # only upper triangle to avoid duplicate comparisons

# Filter out Inf and -Inf values
plot_data_B <- plot_data_B %>%
  filter(is.finite(log2_or))  # Remove rows with Inf or -Inf

# Plot for facet by grade with matrix-style layout (genotype pair matrix within each grade)
ggplot(plot_data_B, aes(x = geno2, y = geno1, fill = log2_or)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", log2_or)), color = "white", size = 4) +
  facet_wrap(~grade, scales = "free") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  labs(
    title = "Plot B: Genotype Pair Differences by Grade (Matrix Layout)",
    x = "Genotype (geno2)",
    y = "Genotype (geno1)",
    fill = "log2(Odds Ratio)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12)
  )












# Prepare data for Plot B (facet by grade, showing all genotype comparisons in a matrix)
plot_data_B <- results_fisher %>%
  filter(!is.na(odds_ratio)) %>%
  mutate(
    log2_or = log2(odds_ratio),
    sig_level = cut(p_adj,
                    breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
                    labels = c("***", "**", "*", "ns")),
    geno_b = factor(geno1, levels = sort(unique(df_props1$resultant_geno))),
    geno_a = factor(geno2, levels = rev(sort(unique(df_props1$resultant_geno)))),
    pair = paste(geno1, "vs", geno2)  # Create pair variable for plotting
  )





# 
# # %>%
#   filter(as.integer(geno1) > as.integer(geno2))  # only upper triangle to avoid duplicate comparisons
# 


# 
# 
# 
# 
# # Prepare data for Plot B (facet by grade, showing all genotype comparisons in a matrix)
# plot_data_B <- results_fisher %>%
#   filter(!is.na(odds_ratio)) %>%
#   mutate(
#     log2_or = log2(odds_ratio),
#     sig_level = cut(p_adj,
#                     breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
#                     labels = c("***", "**", "*", "ns")),
#     geno1 = factor(geno1, levels = rev(sort(unique(df_props1$resultant_geno)))),  # Reversed order for geno1
#     geno2 = factor(geno2, levels = rev(sort(unique(df_props1$resultant_geno)))),  # Reversed order for geno2
#     pair = paste(geno1, "vs", geno2)  # Create pair variable for plotting
#   ) %>%
#   filter(as.integer(geno1) < as.integer(geno2))  # only upper triangle to avoid duplicate comparisons

# Filter out Inf and -Inf values
# plot_data_B <- plot_data_B %>%
#   filter(is.finite(log2_or))  # Remove rows with Inf or -Inf

# Plot for facet by grade with matrix-style layout (genotype pair matrix within each grade)
ggplot(plot_data_B, aes(x = geno_a, y = geno_b, fill = log2_or)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", log2_or)), color = "white", size = 4) +
  facet_wrap(~grade, scales = "free") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  labs(
    title = "Plot B: Genotype Pair Differences by Grade (Matrix Layout)",
    x = "Genotype (geno2)",
    y = "Genotype (geno1)",
    fill = "log2(Odds Ratio)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12)
  )







# Prepare data for Plot B (facet by grade, showing all genotype comparisons in a matrix)
plot_data_B <- results_fisher %>%
  filter(!is.na(odds_ratio)) %>%
  mutate(
    log2_or = log2(odds_ratio),
    sig_level = cut(p_adj,
                    breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
                    labels = c("***", "**", "*", "ns")),
    geno_b = factor(geno1, levels = sort(unique(df_props1$resultant_geno))),
    geno_a = factor(geno2, levels = rev(sort(unique(df_props1$resultant_geno)))),
    pair = paste(geno1, "vs", geno2)  # Create pair variable for plotting
  )



plot_data_B <- plot_data_B %>%
  mutate(sig_level = ifelse(sig_level == "ns", NA, as.character(sig_level)))


# Assuming you have the 'sig_level' column creaNA_character_# Assuming you have the 'sig_level' column created earlier in the 'plot_data_B'
ggplot(plot_data_B, aes(x = geno_a, y = geno_b, fill = log2_or)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", log2_or)), color = "white", size = 4) +
  # Add significance labels based on the 'sig_level' column
  geom_text(aes(label = sig_level), color = "black", size = 5, vjust = -1) +
  facet_wrap(~grade, scales = "free") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  labs(
    title = "Plot B: Genotype Pair Differences by Grade (Matrix Layout)",
    x = "Genotype (geno2)",
    y = "Genotype (geno1)",
    fill = "log2(Odds Ratio)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12)
  )





# Prepare data for Plot B (facet by grade, showing all genotype comparisons in a matrix)
plot_data_B <- results_fisher %>%
  filter(!is.na(odds_ratio)) %>%
  mutate(
    log2_or = log2(odds_ratio),
    sig_level = cut(p_adj,
                    breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
                    labels = c("***", "**", "*", "ns")),
    geno_b = factor(geno1, levels = sort(unique(df_props1$resultant_geno))),
    geno_a = factor(geno2, levels = rev(sort(unique(df_props1$resultant_geno)))),
    pair = paste(geno1, "vs", geno2)  # Create pair variable for plotting
  )
plot_data_B_clean <- plot_data_B %>%
  filter(!is.na(log2_or)) %>%
  mutate(
    sig_level = ifelse(sig_level == "ns", "", as.character(sig_level)),
    label_combined = paste0(sprintf("%.2f", log2_or), "\n", sig_level),
    text_color = ifelse(abs(log2_or) < 3.5, "black", "white")  # adjust threshold if needed
  )

ggplot(plot_data_B_clean, aes(x = geno_a, y = geno_b, fill = log2_or)) +
  geom_tile(color = "black") +  # add black border
  geom_text(
    aes(label = label_combined, color = text_color),
    size = 4.5,
    lineheight = 0.9,
    show.legend = FALSE
  ) +
  scale_color_identity() +
  facet_wrap(~grade, scales = "free") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  labs(
    title = "Plot B: Genotype Pair Differences by Grade (Matrix Layout)",
    x = "Genotype (geno2)",
    y = "Genotype (geno1)",
    fill = "log2(Odds Ratio)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12)
  )
















