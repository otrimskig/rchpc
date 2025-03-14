source("libs.R")

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(gprofiler2)




genes<-readRDS("nf1g/ds/gene_stats.rds")

#all data. 
# gostres<- gost(query = genes$gene_id_ms,
#                organism = "mmusculus",
#                numeric_ns = "ENTREZGENE_ACC")
# 
# 
# gostplot(gostres, capped = F, interactive = F)
# 

#genericise for m-dexps. 


dexp_files<-list.files(path = "m-dexps", pattern= ".rds$", full.names = TRUE)


dexp_files<-dexp_files[1]



#sub("\\.rds$", ".pdf", dexp_files[1])

for (i in 1:length(dexp_files)){
  
i<-1


#html_out_path<-paste0("plots/grof/", "gprof-", basename(sub("\\.rds$", ".html", dexp_files[i])))

d<-readRDS(dexp_files[i])%>%
  filter(FDR<0.01)%>%
  filter(abs(logFC)>2)%>%
  arrange(logFC)


d_name<-basename(dexp_files[i])

}



for (i in 1:length(dexp_files)) {
  
  i <- 1  # Ensure i starts at 1 for debugging
  
  # Load and filter data
  d <- readRDS(dexp_files[i]) %>%
    filter(FDR < 0.01) %>%
    filter(abs(logFC) > 2) %>%
    arrange(logFC)
  
  d_name <- basename(dexp_files[i])
  
  # Run gProfiler analysis
  g <- gost(query = d$gene_id_ms,
            organism = "mmusculus",
            numeric_ns = "ENTREZGENE_ACC",
            ordered_query = TRUE)
  
  # Extract results
  results <- g$result  # Extract enrichment results
 
  significant_labels <- results %>% filter(-log10(p_value) >5)  # Only label p < 0.05
  
  
   
  # Custom ggplot for visualization
  p <- ggplot(results, aes(x = source, y = -log10(p_value), label = term_name)) +
    
    
    geom_point(aes(size = intersection_size, color = source), alpha = 0.8) +  # Color by category (source)
    
  
    
    geom_label_repel(data = significant_labels,  # Only label significant terms
                     aes(fill = source), 
                     size = 3, max.overlaps = 10, show.legend = FALSE) +
    
    
    
    
    # Label pathways
    scale_color_manual(values = c("GO:MF" = "#E41A1C", "GO:CC" = "#377EB8", "GO:BP" = "#4DAF4A", 
                                  "KEGG" = "#FF7F00", "REAC" = "#984EA3", "TF" = "#FFFF33", 
                                  "MIRNA" = "#F781BF", "HP" = "#808080")) +  # Custom colors
    scale_fill_manual(values = c("GO:MF" = "#E41A1C", "GO:CC" = "#377EB8", "GO:BP" = "#4DAF4A", 
                                 "KEGG" = "#FF7F00", "REAC" = "#984EA3", "TF" = "#FFFF33", 
                                 "MIRNA" = "#F781BF", "HP" = "#808080")) +
    labs(title = paste("GProfiler Enrichment:", d_name),
         x = "Term Size", 
         y = "-log10(p-value)",
         color = "Pathway Source") +
    theme_minimal() +
    theme(legend.position = "right")
  
  # Display plot
  print(p)
}

results <- g$result %>%
  mutate(g_label = if_else(-log10(p_value) > 5, 
                           stringr::str_wrap(term_name, width=25),
                           NA))

# Filter only significant terms for labeling
significant_labels <- results %>% filter(!is.na(g_label))


my_pal <- function(range = c(2, 12)) {
  force(range)
  function(x) scales::rescale(x, to = range, from = c(0, 1))
}




p <- ggplot(results, aes(x = source, y = -log10(p_value), label = g_label)) +
 

  

  # Add jittered points
  geom_jitter(aes(size = intersection_size, fill = source, group=source),
              color = "black",
              position = position_jitter(width = 0.2, seed = 3),
              alpha = 0.6, 
              shape = 21,      # Ensures both fill & stroke work
              stroke = 0.5) +  # Controls outline thickness
  
  
  coord_cartesian(clip = "off")

  



plot_data <- ggplot_build(p)$data[[1]]

ggplot(plot_data, aes(x=x, y=y, color=group, size=size, label=label))+
  geom_point()+
  geom_label_repel(aes(label=label))



# Ensure labels match jittered points by setting the same width
  geom_label_repel(aes(label = g_label,
                       point.size = intersection_size), 
                   #xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),
                   point.padding = 0,
                   min.segment.length = 0,
                   size = 3, 
                   max.overlaps = 10, 
                   force = 20,
                   hjust= 0,
                   box.padding = 0,
                   position = position_jitter(width = 0.2, seed=3)) +   
  
  
  
  # geom_text_repel(
  #   data          = subset(dat, wt < 3),
  #   nudge_x       = 2.7 - subset(dat, wt < 3)$wt,
  #   segment.size  = 0.2,
  #   segment.color = "grey50",
  #   direction     = "y",
  #   hjust         = 1
  # ) +  

  


  
  # Scale point sizes with a larger range
  continuous_scale(
    aesthetics = c("size", "point.size"),  # Scale both point size and label size
    palette = my_pal(c(2.5, 12.5))  # Custom palette for scaling
    #guide = guide_legend(override.aes = list(label = ""))  # Hide "a" in legend
  ) +

  # Custom colors for pathway sources
  scale_fill_manual(values = c("GO:MF" = "#E41A1C", "GO:CC" = "#377EB8", "GO:BP" = "#4DAF4A", 
                                "KEGG" = "#FF7F00", "REAC" = "#984EA3", "TF" = "#FFFF33", 
                                "MIRNA" = "#F781BF", "HP" = "#808080")) +  
  
  labs(title = paste("GProfiler Enrichment:", d_name),
       x = "Pathway Source", 
       y = "-log10(p-value)",
       color = "Pathway Source") +
  
  theme_minimal() +
  theme(legend.position = "right")

# Display plot
print(p)

