##########################
source("libs.R")
library(tidyverse)
library(edgeR)
library(dtplyr)
library(ggplot2)
library(ComplexHeatmap)

vm0<-readRDS("acral_sub_rppa/ds/rppa_df_list0.rds")
sample_info0<-readRDS("acral_sub_rppa/ds/sample_info0.rds")
sample_interest<-sample_info0%>%
  filter(sample_type!="none")

vm1<-vm0[["L4 (linear)"]][["df_long"]]
vm2<-vm1%>%
  filter(order %in% sample_interest$order)%>%
  mutate(rppa_value=as.numeric(rppa_value))%>%
  
  left_join(sample_info0%>%select(sample_num, order, sample_type))
  
  # filter(order!=455&order!=456)%>%
  # filter(order!=469)%>%
  
 
saveRDS(vm2, "acral_sub_rppa/ds/all_data_long.rds")
  

# library
library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)
  # 
  # left_join(sample_interest%>%select(mouse_num, order))%>%
  # mutate(sample_name_hm=paste0(mouse_num, "-", order))%>%
  # select(-order, -mouse_num)%>%
  # pivot_wider(names_from = antibody_name, values_from=rppa_value)%>%
  # column_to_rownames("sample_name_hm")%>%
  # as.matrix.data.frame()%>%
  # t()
  # 



rep_diff <- vm2 %>%
  group_by(antibody_name, sample_num, sample_type) %>%
  summarise(
    val1 = rppa_value[1],
    val2 = rppa_value[2],
    mean=mean(c(val1, val2)),
    abs_diff = abs(val1 - val2),
    rel_diff = abs(val1 - val2) / mean(c(val1, val2)),
    .groups = "drop"
  )


vm3<-vm2%>%
  left_join(rep_diff)




# make sure sample_num is a factor for nice ordering
vm3$sample_num <- as.factor(vm3$sample_num)

# plot
p_rel_diff<-ggplot(vm3, aes(x = rel_diff, y = sample_num, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01, fill="black", color="white") +
  #scale_fill_viridis(name = "Relative Diff", option = "E") +
  labs(
    title = "Distribution of Relative Differences by Sample",
    x = "Relative Difference",
    y = "Sample"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )+
  geom_vline(xintercept = 0.375, linetype = "dashed", color = "red", linewidth = 1.2)+
  annotate("text",
           x = 0.375, y = Inf, vjust = 3, hjust=-.05,
           label = "relative difference cutoff", angle = 0,
           size = 7, color = "red")


ggsave("acral_sub_rppa/plots/rel_diff_cutoff.pdf",
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       plot=p_rel_diff,
       limitsize = FALSE,
       
       
       height=2,
       width=2,
       scale = 4
       )


ggplot(vm3, aes(x = rel_diff)) +
  geom_density(linewidth=1,alpha=.2) +
  
  theme_minimal()




filtered_diffs<-vm3%>%
  mutate(rppa_filtered=if_else(rel_diff<=.375, rppa_value, NA))%>%
  mutate(mean_filtered=if_else(rel_diff<=.375, mean, NA))



hm_mat<-filtered_diffs%>%
  mutate(sample_num_type=paste0(sample_num, "_", sample_type))%>%
  select(antibody_name, sample_num_type, mean_filtered)%>%
  unique()%>%
  pivot_wider(names_from = sample_num_type, values_from=mean_filtered)%>%
  column_to_rownames("antibody_name")%>%
  as.matrix.data.frame()


saveRDS(hm_mat, "acral_sub_rppa/ds/umat-per_sample_375_rel_diff_th.rds")




n_threshold <- 9  # for example, keep rows with at least 5 non-NA entries

# Apply filtering to the matrix
hm_mat_filtered <- t(scale(t(hm_mat[rowSums(!is.na(hm_mat)) >= n_threshold, ])))



Heatmap(hm_mat_filtered)





fdf2<-filtered_diffs%>%
  select(-rppa_value, -rppa_filtered, -order)%>%
  unique()
  
library(dplyr)
library(DESeq2) # For the exact test



fdf2%>%
  arrange(mean)%>%
  str()


library(dplyr)
library(purrr)

# Split by antibody_name
test_results <- fdf2 %>%
  group_by(antibody_name) %>%
  group_split() %>%
  map_dfr(function(df) {
    ab_name <- unique(df$antibody_name)
    groups <- df %>%
      group_by(sample_type) %>%
      summarise(mean_val = mean(mean, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(mean_val))
    
    n_groups <- nrow(groups)
    
    # T-test (2 groups)
    t_test_p <- if (n_groups == 2) {
      tryCatch(t.test(mean ~ sample_type, data = df)$p.value, error = function(e) NA)
    } else NA
    
    # Fold-change (between top 2 sample types by mean)
    fold_change <- if (n_groups >= 2) {
      fc <- groups$mean_val[1] / groups$mean_val[2]
      if (is.finite(fc)) fc else NA
    } else NA
    

    
    tibble(
      antibody_name = ab_name,
      t_test_p = t_test_p,
     
      fold_change = fold_change
    )
  })






test_results <- fdf2 %>%
  group_by(antibody_name) %>%
  group_split() %>%
  map_dfr(function(df) {
    ab_name <- unique(df$antibody_name)
    
    # Summary stats per group
    group_stats <- df %>%
      group_by(sample_type) %>%
      summarise(
        n = n(),
        mean_val = mean(mean, na.rm = TRUE),
        sd_val = sd(mean, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(mean_val))
    
    n_groups <- nrow(group_stats)
    
    # T-test (2 groups)
    t_test_p <- if (n_groups == 2) {
      tryCatch(t.test(mean ~ sample_type, data = df)$p.value, error = function(e) NA)
    } else NA
    
    # Fold-change (between top 2 sample types by mean)
    fold_change <- if (n_groups >= 2) {
      fc <- group_stats$mean_val[1] / group_stats$mean_val[2]
      if (is.finite(fc)) fc else NA
    } else NA
    
    # Spread out stats using sample_type names as prefixes
    group_cols <- group_stats %>%
      pivot_longer(cols = c(n, mean_val, sd_val)) %>%
      mutate(name = paste0(sample_type, "_", name)) %>%
      select(name, value) %>%
      pivot_wider(names_from = name, values_from = value)
    
    tibble(
      antibody_name = ab_name,
      t_test_p = t_test_p,
      fold_change = fold_change
    ) %>%
      bind_cols(group_cols)
  })







#plot t-test distributions.
ggplot(test_results, aes(x = t_test_p)) +
  geom_density(linewidth=1,alpha=.2)





test_results2 <- fdf2 %>%
  group_by(antibody_name) %>%
  group_split() %>%
  map_dfr(function(df) {
    ab_name2 <- unique(df$antibody_name)
    
    # Summary stats per group, using actual sample_type names
    group_stats <- df %>%
      group_by(sample_type) %>%
      summarise(
        n = sum(!is.na(mean_filtered)),
        mean_val = mean(mean_filtered, na.rm = TRUE),
        sd_val = sd(mean_filtered, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(mean_val))
    
    n_groups <- nrow(group_stats)
    
    # T-test (2 groups only)
    t_test_p2 <- if (n_groups == 2) {
      tryCatch(t.test(mean_filtered ~ sample_type, data = df)$p.value, error = function(e) NA)
    } else NA
    
    # Fold-change (between top 2 sample types by mean)
    fold_change2 <- if (n_groups >= 2) {
      fc <- group_stats$mean_val[1] / group_stats$mean_val[2]
      if (is.finite(fc)) fc else NA
    } else NA
    
    # Spread out stats using actual sample_type names as prefixes
    group_cols <- group_stats %>%
      pivot_longer(cols = c(n, mean_val, sd_val)) %>%
      mutate(name = paste0(sample_type, "_", name)) %>%
      select(name, value) %>%
      pivot_wider(names_from = name, values_from = value)
    
    tibble(
      antibody_name = ab_name2,
      t_test_p = t_test_p2,
      fold_change = fold_change2
    ) %>%
      bind_cols(group_cols)
  })



#plot t-test distributions.
ggplot(test_results2, aes(x = t_test_p)) +
  geom_density(linewidth=1,alpha=.2)




test_results2.1 <- test_results2 %>%
  select(sort(names(.))) %>%
  relocate(antibody_name, t_test_p) %>%
  mutate(across(everything(), ~ ifelse(is.nan(.), NA, .))) %>%
  rename_with(~ ifelse(.x == "antibody_name", .x, paste0(.x, "_excl")))




saveRDS(test_results2.1, "acral_sub_rppa/ds/stats_fr_threshold_var.rds")
saveRDS(test_results, "acral_sub_rppa/ds/stats_no_exclusion.rds")
