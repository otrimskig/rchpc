source("libs.R")
suppressPackageStartupMessages(library(ComplexHeatmap))
# ht_opt$message = FALSE
library(gridtext)
library(grid)
library(gridExtra)

library(circlize)
library(purrr)
library(tidyverse)
library(ggplot2)
library(crayon)
library(fs)
library(Cairo)


source("ggplot_aspect_function.R")
library(dplyr)
library(stringr)

# Function to identify gene names (assuming gene names are uppercase words)
is_gene_name <- function(x) {
  # A simple pattern to detect gene names, assuming they are uppercase with numbers
  return(str_detect(x, "^[A-Z0-9]+$"))
}

# Function to capitalize pathway names, leaving gene names in uppercase
capitalize_pathway <- function(x) {
  words <- unlist(str_split(x, " "))  # Split the text into words
  
  # Capitalize each word unless it's recognized as a gene name
  words <- sapply(words, function(word) {
    if (is_gene_name(word)) {
      return(word)  # Leave gene names in uppercase
    } else {
      return(str_to_title(word))  # Title-case the rest
    }
  })
  
  # Recombine the words into a single string
  return(paste(words, collapse = " "))
}


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
  dplyr::filter(!grepl("MALE", Pathway))%>%
  
  mutate(Pathway_full=Pathway)%>%
  mutate(Pathway = stringr::str_trunc(Pathway, width = 50))%>%
  mutate(Pathway = make.unique(Pathway))






df0<-dfu1%>%
  pivot_longer(cols= -c(Pathway, Pathway_full), names_to = "mouse_num", values_to = "z_score")%>%
  left_join(sample_info%>%
              select(mouse_num, patho_cat_name, resultant_geno))%>%
  
  left_join(analysis_list$results2)%>%
  
  mutate(Pathway=gsub("_", " ", Pathway))%>%
  
  mutate(Pathway=sub("h\\.\\.", "h\\.", Pathway))%>%
  mutate(Pathway = str_to_title(Pathway))%>%

  
  mutate(Pathway = as.character(Pathway))

df1 <- df0 %>%
  filter(
    min_pval <= quantile(min_pval, 0.05, na.rm = TRUE),  # lowest 10% p-values
    abs_diff >= quantile(abs_diff, 0.95, na.rm = TRUE)   # top 10% abs_diff
  )%>%
  filter(min_pval<=.05)



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
  
  
  # Step 1: Original and padded labels
  original_labels <- pathway_chunks[[i]]
  padded_labels <- pad_labels(original_labels, width = 50)
  
  # Step 2: Create named map and preserve order
  label_map <- setNames(padded_labels, original_labels)
  ordered_padded_labels <- unname(label_map[original_labels])  # Preserve the ordering!
  
  # Step 3: Create plot data
  df_chunk <- df2 %>%
    filter(Pathway %in% original_labels) %>%
    mutate(
      Pathway = factor(Pathway, levels = rev(original_labels)),  # controls grouping
      
      # Display label with fixed width formatting
      Pathway_label = factor(label_map[as.character(Pathway)],
                             levels = rev(ordered_padded_labels))  # controls label order
    )

  
  p <- ggplot(data = df_chunk, aes(x = z_score, y = Pathway_label, color = resultant_geno)) + 
    
    geom_segment(aes(x = mean_group_a, xend = mean_group_b,
                     y = Pathway_label, yend = Pathway_label),
                 color="gray", alpha=.1, linewidth=2)+
    
    geom_point(aes(x = mean_group_a, y = Pathway_label), color = "#B79F00", shape = 124, size = 3) +
    geom_point(aes(x = mean_group_b, y = Pathway_label), color = "#00BFC4", shape = 124, size = 3) +
    geom_point(size = 2, alpha = 0.35) +
   
   
    ggh4x::facet_nested(pathway_grouping ~ patho_cat_name, scales = "free_y", space = "free") +
    theme_bw() +
    theme(
      strip.placement = "outside",
      strip.text = element_text(size = 10),
      strip.background = element_blank(),

      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      panel.spacing = unit(0.0, "lines"),
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(100, 100, 100, 100),
      axis.text.y = element_text(size = 7) #family = "Arial", 
      #legend.position = "right",         # (x, y) coordinates in NPC units (0–1)
      #legend.justification = "top",    # anchor the legend box to its top-right corner
      #legend.direction = "vertical"     # or "horizontal" depending on layout
    ) +
    
    labs(color = "Resultant Genotype")+
    
    scale_x_continuous(limits = c(-.5,.5))+

    scale_y_discrete(expand = expansion(add = c(1, 1))) +
    scale_color_manual(values = col_map[["resultant_geno"]])
    
  
  p2<-p
  
 # determine plot aspect ratio dynamically to make it consistent

  asp<-gg_output_aspect(
    p,
    panel_width = 1,
    row_height = .3,
    n_rows = length(pathway_chunks[[i]]) # Number of unique y-axis rows (across facets)
  )

  # CairoPDF(
  #   file = paste0("nf1g/gsvas/plots/chunks_per_patho/", 
  #                 name_of_patho_cat, 
  #                 sprintf("plot_chunk_%02d.pdf", i)),
  #   width = 1.2,
  #   height = 1 / asp,
  #   pointsize = 12  # Adjust point size if necessary
  # )
  # print(p)  # Plot your ggplot object
  # dev.off()  # Close the device
 
  
  # # Save plot
  ggsave(
    filename = paste0("nf1g/gsvas/plots/chunks_per_patho/",
                       name_of_patho_cat,
                       sprintf("plot_chunk_%02d.pdf", i)),
    device="pdf",
    plot = p2,
    width = 1.2,
    dpi=600,
    #height = 10,
    height = 1/asp,
    scale = 10,
    limitsize = FALSE
  )


  
  
  
  cat(cli::col_blue(i), "of", cli::col_red(length(pathway_chunks)), "\n")
}





}





