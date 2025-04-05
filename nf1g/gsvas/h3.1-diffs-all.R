source("libs.R")
# suppressPackageStartupMessages(library(ComplexHeatmap))
# ht_opt$message = FALSE
library(gridtext)
library(circlize)
library(purrr)
library(tidyverse)
library(ggplot2)
library(crayon)
library(fs)


#get all info to be used across loops.
sample_info<-readRDS("nf1g/ds/v10-per_sample_updated.rds")
matu0<-readRDS("timex/ds/hu-msig_all-gsva-values.rds")
col_map<-readRDS("nf1g/surv/colors_map_surv.rds")



#get paths of rds file analyses.
a<-dir_ls("nf1g/gsvas/ds")%>%
  tibble(path=.)%>%
  filter(grepl("a86", path))%>%
  filter(grepl("nf1 KO; pten KO; ink KO; atrx KOvsnf1 KO; pten KO; ink KO; atrx w", path))

rds_diff_files<-a$path




#start loops here.

pa<-1



for(pa in 1:length(rds_diff_files)){


analysis_list<-readRDS(rds_diff_files[pa])

name_of_patho_cat<-sub(" :.*","",analysis_list$metadata$group_comparison)


dfu1<-matu0[,c(analysis_list$metadata$group_a_samples, 
                analysis_list$metadata$group_b_samples), drop=FALSE]%>%
  as_tibble(rownames="Pathway")%>%

  dplyr::filter(!grepl("^c1", Pathway))%>%
  dplyr::filter(!grepl("^c3.MIR", Pathway))%>%
  dplyr::filter(!grepl("SEX", Pathway))%>%
  dplyr::filter(!grepl("GENDER", Pathway))%>%
  dplyr::filter(!grepl("MALE", Pathway))






df0<-dfu1%>%
  pivot_longer(cols= -Pathway, names_to = "mouse_num", values_to = "z_score")%>%
  left_join(sample_info%>%
              select(mouse_num, patho_cat_name, resultant_geno))%>%
  
  left_join(analysis_list$results2)%>%
  
  mutate(Pathway=gsub("_", " ", Pathway))%>%
  mutate(Pathway=gsub("\\.\\.", "\\.", Pathway))%>%
  
  mutate(Pathway = as.character(Pathway))

df1 <- df0 %>%
  filter(
    min_pval <= quantile(min_pval, 0.05, na.rm = TRUE),  # lowest 10% p-values
    abs_diff >= quantile(abs_diff, 0.95, na.rm = TRUE)   # top 10% abs_diff
  )



df2<-df1%>%
  mutate(pathway_grouping=substr(Pathway, 1,2))%>%
  mutate(pathway_grouping=sub("^h\\.$", "h", pathway_grouping))%>%
  select(-contains("No evidence"))%>%
  filter(!grepl("^c1", Pathway))%>%
  filter(!grepl("^c3.MIR", Pathway))%>%
  mutate(diff_direction=sign(diff))

# Step 1: Order by pathway_grouping and descending diffs
ordered_pathways <- df2 %>%
  distinct(Pathway, pathway_grouping, diff, mean_total, diff_direction) %>%
  arrange(pathway_grouping, desc(diff_direction), desc(mean_total)) %>%
  pull(Pathway)




# Step 2: Create chunks of 1000 based on the above order
pathway_chunks <- split(ordered_pathways, ceiling(seq_along(ordered_pathways) / 1000))



# Step 3: Output directory
dir.create("nf1g/gsvas/plots/chunks_per_patho", recursive = TRUE, showWarnings = FALSE)






# Loop through and plot each chunk
for (i in seq_along(pathway_chunks)) {
  
  df_chunk <- df2 %>%
    filter(Pathway %in% pathway_chunks[[i]]) %>%
    mutate(Pathway = factor(Pathway, levels = rev(pathway_chunks[[i]])))  # reversed order for top-down
  
  p <- ggplot(data = df_chunk, aes(x = z_score, y = Pathway, color = resultant_geno)) + 
    geom_point(aes(x = mean_group_a, y = Pathway), color = "black", shape = 124, size = 3) +
    geom_point(aes(x = mean_group_b, y = Pathway), color = "black", shape = 124, size = 3) +
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
    scale_color_manual(values = col_map[["resultant_geno"]])
  
  # Save plot
  ggsave(
    filename = paste0("nf1g/gsvas/plots/chunks_per_patho/", 
                       name_of_patho_cat,
                       sprintf("plot_chunk_%02d.pdf", i)),
    plot = p,
    width = 2,
    height = 4,
    scale = 8,
    limitsize = FALSE
  )
  
  cat(cli::col_blue(i), "of", cli::col_red(length(pathway_chunks)), "\n")
}





}











