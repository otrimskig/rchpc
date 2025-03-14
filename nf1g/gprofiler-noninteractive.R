source("libs.R")

library(tidyverse)
library(gprofiler2)




genes<-readRDS("nf1g/ds/gene_stats.rds")

#all data. 
gostres<- gost(query = genes$gene_id_ms,
               organism = "mmusculus",
               numeric_ns = "ENTREZGENE_ACC")


gostplot(gostres, capped = F, interactive = F)


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



g<- gost(query = d$gene_id_ms,
               organism = "mmusculus",
               numeric_ns = "ENTREZGENE_ACC",
         ordered_query = TRUE)

p<-gostplot(g, capped = TRUE, interactive = F)




#save as self-contained html file. 
# htmlwidgets::saveWidget(p, html_out_path)



}



p

# publish_gosttable(g)




library(ggplot2)
library(ggrepel)
library(dplyr)
library(gprofiler2)

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
  mutate(g_label = if_else(-log10(p_value) > 5, term_name, NA))

# Filter only significant terms for labeling
significant_labels <- results %>% filter(!is.na(g_label))

p <- ggplot(results, aes(x = source, y = -log10(p_value))) +
 
  
  
  
  
  
  
  
  
  # Scale point sizes with a larger range
  scale_size_continuous(range = c(2, 12)) +  # Adjust range as needed
  
  
  
  
  
  
  # Ensure labels match jittered points by setting the same width
  geom_label_repel(aes(label = g_label), 
                   size = 3, 
                   max.overlaps = 10, 
                   force = 50,
                   hjust= 0,
       
                   position = position_jitter(width = 0.2, seed=3)) + 
  
  
  
  
  # Add jittered points
  geom_jitter(aes(size = intersection_size), 
              fill = "white",
              color = "black",
              position = position_jitter(width = 0.2, seed = 3),
              alpha = 1, 
              shape = 21,      # Ensures both fill & stroke work
              stroke = 0) +  # Controls outline thickness
  
  # Add jittered points
  geom_jitter(aes(size = intersection_size, fill = source),
              color = "black",
              position = position_jitter(width = 0.2, seed = 3),
              alpha = 0.6, 
              shape = 21,      # Ensures both fill & stroke work
              stroke = 0.5) +  # Controls outline thickness
  
  
  
  
  
  




  
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

