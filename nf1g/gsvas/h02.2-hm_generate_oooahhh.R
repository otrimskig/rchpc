source("libs.R")
# suppressPackageStartupMessages(library(ComplexHeatmap))
# ht_opt$message = FALSE
library(gridtext)
library(circlize)
library(purrr)
library(tidyverse)
library(ggplot2)
library(crayon)

sample_info<-readRDS("nf1g/ds/v10-per_sample_updated.rds")
matu0<-readRDS("timex/ds/hu-msig_all-gsva-values.rds")
col_map<-readRDS("nf1g/surv/colors_map_surv.rds")

stats_files<-list.files("nf1g/gsvas/ds", pattern="all_genos", full.names=T)

all_info<-list()
for (fi in 1:length(stats_files)){
  
  all_info[[fi]]<-readRDS(stats_files[fi])
  
}


mean_tib<-tibble(Pathway=NA)

for(em in 1:length(all_info)){

# Extract group_a_name from metadata
group_a_name <- all_info[[em]]$metadata$group_a
group_b_name <- all_info[[em]]$metadata$group_b

# Select Pathway and mean_group_a columns, then rename mean_group_a
mean_tib<-all_info[[em]]$results %>%
  select(Pathway, mean_group_a, mean_group_b,
         std_group_a, std_group_b) %>%
  rename(!!paste0(group_a_name, "_mean") := mean_group_a)%>%
  rename(!!paste0(group_b_name, "_mean") := mean_group_b)%>%
  
  rename(!!paste0(group_a_name, "_sd") := std_group_a)%>%
  rename(!!paste0(group_b_name, "_sd") := std_group_b)%>%
  
  full_join(mean_tib)

}


saveRDS(mean_tib, "nf1g/gsvas/ds/hu-gsva_means_across_genos.rds")


mean_tib2<-mean_tib%>%
  filter(!is.na(Pathway))%>%

  mutate(disease_mean = rowMeans(select(., `Spindle and epithelioid tumors_mean`, 
                                        `Nerve sheath tumors_mean`, 
                                        `Glioneuronal tumors_mean`, 
                                        Gliomas_mean), na.rm = TRUE))%>%
  mutate(abs_disease_mean=abs(disease_mean))%>%
  mutate(mean_diff=disease_mean-`No evidence of disease_mean`)%>%
  mutate(abs_mean_diff=abs(mean_diff))%>%
  
  mutate(pathway_grouping=substr(Pathway, 1,2))%>%
  mutate(pathway_grouping=sub("^h\\.$", "h", pathway_grouping))%>%
  select(-contains("No evidence"))%>%
filter(!grepl("^c1", Pathway))%>%
filter(!grepl("^c3.MIR", Pathway))

  
  
  # 
  # 
  # 

  
  
  
  
  # filter(grepl("NEURO", Pathway)|
  #          grepl("GLIO", Pathway))
 


n <- 2000000  # Choose the number of top pathways per column

# Identify numeric columns to rank
numeric_cols <- mean_tib2 %>%
  select(where(is.numeric)) %>%
  select(contains("_mean"))%>%
  select(-starts_with("abs_"))%>%
  colnames()

# # Get the top `n` pathways for each column
# top_indices <- numeric_cols %>%
#   lapply(function(col) mean_tib2 %>%
#            arrange(desc(.data[[col]])) %>%
#            slice_head(n = n) %>%
#            pull(Pathway)) %>%
#   unlist() %>%
#   unique()  # Ensure unique pathways
# 
# # Subset the main tibble
# mean_tib_top <- mean_tib2 %>%
#   filter(Pathway %in% top_indices)
# 
# 


top_indices <- numeric_cols %>%
  lapply(function(col) {
    mean_tib2 %>%
      group_by(pathway_grouping) %>%
      arrange(desc(.data[[col]]), .by_group = TRUE) %>%
      slice_head(n = n) %>%
      pull(Pathway)
  }) %>%
  unlist() %>%
  unique()

# Subset to keep only the top pathways
mean_tib_top <- mean_tib2 %>%
  filter(Pathway %in% top_indices)



df1<-matu0[mean_tib_top$Pathway,, drop=FALSE]%>%
  as_tibble(rownames="Pathway")%>%
  pivot_longer(cols= -Pathway, names_to = "mouse_num", values_to = "z_score")%>%
  left_join(sample_info%>%
              select(mouse_num, patho_cat_name, resultant_geno))%>%
  left_join(mean_tib_top)%>%
  filter(patho_cat_name!="No evidence of disease")




df_means <- df1 %>%
  group_by(patho_cat_name, Pathway) %>%
  summarize(mean_z_score = mean(z_score, na.rm = TRUE))

df2<-df1%>%
  left_join(df_means)%>%
  mutate(Pathway=gsub("_", " ", Pathway))%>%
  mutate(Pathway=gsub("\\.\\.", "\\.", Pathway))


# ggplot(data = df2, aes(x = z_score, y =reorder(Pathway, disease_mean), color = patho_cat_name)) + 
#   
#   
#   # Add mean points for each group (using black color and vertical bar shape)
#   geom_point(aes(x = mean_z_score, y = Pathway), color = "black", shape = 124, size = 4) +
#   geom_point(size = 2.5, alpha = 0.5) +
#   
#   facet_grid(~pathway_grouping~patho_cat_name, scales = "free_y") +  # Facet by patho_cat_name
#   theme_bw() +
#   theme(
#     strip.placement = "outside",
#     strip.text = element_text(size = 10),
#     strip.background = element_blank(),
#     axis.text.y = element_text(size = 5),  # Larger y-axis labels
#     axis.ticks.y = element_blank(),
#     axis.title.y = element_blank(),
#     panel.spacing = unit(0, "lines"),
#     panel.grid.major.y = element_line(color = "grey90"),
#     panel.grid.minor.y = element_blank(),
#     plot.margin = margin(10, 10, 10, 10)  # Add more padding around the plot
#   ) +
#   scale_y_discrete(expand = expansion(add = c(1, 1))) +
# 
#   #guides(color = "none")
# scale_color_manual(values = col_map[["hist_cat_name"]])



library(ggrastr)
library(ggh4x)

p<-ggplot(data = df2, aes(x = z_score, y = reorder(Pathway, disease_mean), color = patho_cat_name)) + 
  
  geom_point(aes(x = mean_z_score, y = Pathway), 
                       color = "black", shape = 124, size = 3) +
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
  scale_color_manual(values = col_map[["hist_cat_name"]])


plot_object<-p

ggsave("nf1g/gsvas/plots/overall_plot.pdf",
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       plot=plot_object,
       limitsize = FALSE,
       
       
       height=50,
       width=.2,
       scale = 100,
       dpi=10,
       
       
       
)









##################


# 
# # Subset the matrix using these indices
# top_diff <- comp_mat%>%
#   as_tibble()%>%
#   #filter(min_pval<.05)%>%
#   arrange(desc(abs_diff))%>%
#   slice(1:50)
# 
# 
# 
# 
# diff_df1 <- sa_mat0[top_diff$Pathway, , drop = FALSE] %>%
#   as_tibble(rownames = "Pathway")%>%
#   pivot_longer(cols= -Pathway, names_to = "mouse_num", values_to = "z_score")%>%
#   left_join(sample_info%>%
#               select(mouse_num, patho_cat_name, resultant_geno))
#   
# 
# 
# diff_df2<-diff_df1%>%
#   left_join(top_diff%>%select(Pathway, abs_diff, mean_total, diff))
# 
# ggplot(data=diff_df2, aes(x=z_score, y=reorder(Pathway, diff)))+
#   geom_point(aes(color=resultant_geno), size=3, alpha=.5)+
#   theme_bw()
# 
