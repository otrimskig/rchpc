source("libs.R")
library(tidyverse)
library(dtplyr)
suppressPackageStartupMessages(library(ComplexHeatmap))
ht_opt$message = FALSE
library(gridtext)
library(circlize)
library(purrr)




sample_meta<-readRDS("nf1g/ds/v10-per_sample_updated.rds")%>%
  filter(patho_cat!="NED")
  
  

counts<-readRDS("nf1g/ds/vm-02-filtered_rpkms.rds")



df0<-counts%>%
  semi_join(sample_meta%>%
              select(sample_id))%>%
  left_join(sample_meta)%>%

  select(gene_id_ms, gene_name_ms, mouse_num, rpkm)%>%
  filter(!is.na(gene_name_ms))



df1<-df0%>%
  select(-gene_id_ms)%>%
  pivot_wider(names_from = mouse_num, values_from = rpkm)%>%
  group_by(gene_name_ms)%>%
  slice(1)%>%
  ungroup()


mat0<-df1%>%
  column_to_rownames("gene_name_ms")%>%
  as.matrix.data.frame()


matz0<-t(scale(t(mat0)))


matz1<-matz0[1:1000,]

matz1<-matz0


set.seed(889)
sub <- sample(rownames(matz0), size = 1000)
matz1 <- matz0[sub,
                  ]






de_anno_df<-sample_meta%>%select(patho_grade, patho_cat_name, 
                                #patho_cat2_name,
                               # patho_cat_det_name,
                               resultant_geno)

anno_color_set<-readRDS("nf1g/ds/colors_list.rds")
de_anno_colnames<-colnames(de_anno_df)
anno_subset<-list()

for (colnum in 1:length(de_anno_colnames)){
  
  #set column num(this wil be dynamic in for loop.)
  # colnum<-1
  
  #get column name as character, from dynamically defined number.
  colname<-de_anno_colnames[colnum]
  
  #from that column name, get all the unique values.
  unique_to_filter<-de_anno_df%>%select(!!sym(colname))%>%unique()%>%pull()%>%as.character()
  
  
  
  
  anno_subset[[colname]]<-anno_color_set[[colname]][c(unique_to_filter)]
  
  
}


anno<-HeatmapAnnotation(df=de_anno_df, 
                        col=c(anno_subset),
                        annotation_name_side = "left",
                        
                        gp = gpar(col = "black", fontsize = 3),
                        
                        #annotation_label = c(col_proper_names),
                        annotation_name_gp= gpar(fontsize=7),
                        
                        simple_anno_size = unit(.125, "in"))






Heatmap(matz1,
        top_annotation = anno,
        show_row_names = FALSE,
        raster_by_magick = FALSE
        )



# Your matrix is matz1
col_clust <- hclust(dist(t(matz0)))



col_dend <- as.dendrogram(col_clust)



library(dendextend)

# Create dendrogram
col_dend <- as.dendrogram(hclust(dist(t(matz1))))

# Example annotation: color labels by a factor
# Suppose you have some group labels for each column
col_groups <- sample(c("A", "B"), ncol(matz1), replace = TRUE) # just an example

# Color labels by group
labels_colors(col_dend) <- ifelse(col_groups == "A", "red", "blue")

# Plot
plot(col_dend, main = "Column Clustering Dendrogram")





library(ggdendro)

# Create dendrogram
col_clust <- hclust(dist(t(matz1)))

# Convert to ggdendro object
dendro_data <- dendro_data(col_clust)

# Plot using ggplot2
library(ggplot2)
ggplot(segment(dendro_data)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  theme_minimal() +
  #coord_flip() +           # Flip to have the dendrogram vertical
  #scale_y_reverse() +       # Reverse y axis if you want
  labs(title = "Column Dendrogram")


library(ggplot2)
library(ggdendro)

# 1. Cluster the columns
col_clust <- hclust(dist(t(matz1)))

# 2. Cut into k clusters
k <- 8
clusters <- cutree(col_clust, k = k)

# 3. Get dendrogram data for ggplot
dendro_data <- dendro_data(col_clust, type = "rectangle")  # can be 'rectangle' or 'triangle'

# 4. Create dataframe for tip labels
labels_df <- dendro_data$labels
labels_df$cluster <- clusters[labels_df$label]  # match clusters to labels

# 5. Plot the dendrogram
ggplot() +
  geom_segment(data = dendro_data$segments, 
               aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = labels_df, 
             aes(x = x, y = 0, color = factor(cluster)), 
             size = 4) +  # big points at tips colored by cluster
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  coord_flip() +             # Flip x and y axes for vertical dendrogram
  scale_y_reverse() +        # Reverse y-axis so dendrogram points downward
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()) +
  labs(title = "Column Dendrogram with Cluster Annotations", color = "Cluster")







dend_mat <-t(mat0) %>% dist %>% hclust(method = "average") %>% as.dendrogram



# For each annotation column, map the values to colors
anno_colors_matrix <- sapply(names(anno_subset), function(var) {
  var_values <- de_anno_df[[var]]   # e.g., "Drug", "Control"
  color_map <- anno_subset[[var]]   # named colors
  unname(color_map[var_values])     # match values to colors
})

# Transpose if needed to match colored_bars() expectation
#anno_colors_matrix <- t(anno_colors_matrix)


par(mar = c(5,2,2,1))

labels_colors(col_dend) <- "black"
plot(col_dend, main = "Column Dendrogram with Annotations")

colored_bars(colors = anno_colors_matrix, 
             dend = col_dend, 
             rowLabels = names(anno_subset),
             sort_by_labels_order = TRUE)




ggd1 <- as.ggdend(col_dend)


ggplot(ggd1)
