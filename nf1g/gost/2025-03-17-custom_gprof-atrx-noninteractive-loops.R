source("libs.R")

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(gprofiler2)
library(ggpp)
library(crayon)
library(grid)
library(gridExtra)
library(colorspace)

set.seed(42)
gost_v1<-readRDS("nf1g/atrx_comps/ds/gost_v1.rds")
#gost_v1<-gost_v1[3]



pathway_list_sum<-tibble(pathway_list=NA, index=NA, index_index=NA, p_value=NA)
for (ps in 1:length(gost_v1)){
  
  pathway_list_sum<-bind_rows(pathway_list_sum,
                              as_tibble(gost_v1[[ps]]$result)%>%
                                mutate(index=1:n())%>%
                                mutate(index_index=ps)%>%
                                rename(pathway_list=source)%>%
                                dplyr::select(pathway_list, index, index_index, p_value))%>%
    drop_na()

}





pathway_list_tib<-tibble(pathway_list=NA)
for (de1.5 in 1:length(gost_v1)){
  
  pathway_list_tib<-full_join(pathway_list_tib,
                              as_tibble(gost_v1[[de1.5]]$result)%>%
                                rename(pathway_list=source)%>%
                                dplyr::select(pathway_list)%>%unique())%>%
    unique()%>%
    drop_na()%>%
    arrange(pathway_list)
  
  pathway_list_factor<-pathway_list_tib%>%
    mutate(pathway_list=factor(pathway_list, levels=pathway_list_tib$pathway_list))%>%
    mutate(group_x_coord=1:n())
  
  
}

# 
# 
# pathway_palette <- setNames(
#   scales::hue_pal()(length(pathway_list_factor$pathway_list)),# Ensure this matches pathway_list categories
#   pathway_list_factor$pathway_list
# )
# 
# 
# 
# # Function to lighten color by blending with gray
# lighten_color <- function(col, factor = 0.5) {
#   col_rgb <- col2rgb(col) / 255  # Normalize RGB values
#   gray_rgb <- col_rgb * (1 - factor) + factor  # Blend with gray
#   rgb(gray_rgb[1], gray_rgb[2], gray_rgb[3])  # Convert back to hex
# }
# 
# # Create a modified version of the palette where high p-values are lightened
# adjusted_palette <- pathway_palette  # Copy original colors
# 
# 











for (de2 in 1:length(gost_v1)){

# de2<-3
  
cat(blue("working on "), red(names(gost_v1$de2)), " (", cyan(de2), " of ", length(gost_v1), ")", "\n")
  

data1<-gost_v1[[de2]]$result

data2 <- data1%>%
  mutate(point_index=1:n())%>%
  rename(pathway_list=source)

centered_jitter <- function(x, width = 0.65, seed = NULL) {
  if (!is.null(seed)) set.seed(seed) # Ensure reproducibility
  n <- length(x)
  if (n == 1) {
    jitter_offsets<-0
    return(x + jitter_offsets)  # Keep single points unchanged
  } else {
    jitter_offsets <- runif(n, min = -width / 2, max = width / 2)
    jitter_offsets <- jitter_offsets - mean(jitter_offsets)  # Center around 0
    return(x + jitter_offsets)
  }
}


# Apply to dataset
data3 <- left_join(pathway_list_factor,data2)%>%
  
  group_by(pathway_list)%>%
  mutate(jittered_x = centered_jitter(group_x_coord, width = 0.4, seed = 3))%>%
  mutate(dist_from_group_x=jittered_x-group_x_coord)%>%
  
  mutate(label_x_pos=if_else(dist_from_group_x<0, group_x_coord-.3, group_x_coord+.3))%>%
  
  group_by(pathway_list) %>% #create ranked list for filtering labels.
  mutate(
    group_rank = rank(p_value, ties.method = "average"),  # Rank p-values within each pathway
    group_percentile = (1 - percent_rank(p_value)) * 100  # Convert to percentile (0-100 scale)
  ) %>%
  ungroup()%>%
  mutate(
    rank = rank(p_value, ties.method = "average"),  # Rank p-values within each pathway
    percentile = (1 - percent_rank(p_value)) * 100  # Convert to percentile (0-100 scale)
  )%>%
  

######set labelling threshold############################################  
  mutate(g_label = if_else(-log10(p_value)>5 &
                            ((grepl("MHC", term_name)&group_percentile>=80) |group_rank<=5|group_percentile>=97),
                              
                           stringr::str_wrap(term_name, width=25),
                           
                           NA))


########################################################



# 
# high_pvals <- data3$p_value > 5   # Define high p-values
# 
# # Apply lightening only to pathways with high p-values
# adjusted_palette[data3$pathway_list[high_pvals]] <- 
#   sapply(pathway_palette[data3$pathway_list[high_pvals]], lighten_color, factor = 0.5)
# 
# 
# 



my_pal <- function(range = c(1, 10)) {
  force(range)
  function(x) scales::rescale(x, to = range, from = c(0, 1))
}

x_scale_factor_man<-.5
breaks_dynamic<-seq(x_scale_factor_man, by=x_scale_factor_man, length.out=length(pathway_list_factor$pathway_list))
limits_dynamic<-c(min(breaks_dynamic)-.25,
                  max(breaks_dynamic)+.25)




set.seed(69)
gost_v1[[de2]]$plot1 <- ggplot(data3, aes(x = group_x_coord*x_scale_factor_man+dist_from_group_x, 
                        y = -log10(p_value), 
                        label = g_label)) +
  
  # Add jittered points

  
  # Add jittered points
  geom_point(aes(size = intersection_size, 
                 color = pathway_list # Lighten color
                 ),
             shape = 16,
             alpha = 0.4)+# Ensures both fill & stroke work) +      
  geom_point(aes(size = intersection_size, 
                 color = pathway_list # Lighten color
                 ),
             shape = 1,
             stroke=.8,
             alpha = 0.8)+# Ensures both fill & stroke work) +  
  
  
  # scale_color_manual(values = setNames(
  #   desaturate(unique(data3$pathway_list), 0.5),  # Precompute desaturated colors
  #   unique(data3$pathway_list)
  # )) +

  geom_label_repel(data=subset(data3),
                   aes(label = g_label,
                       point.size = intersection_size,
                       fill=pathway_list),
                   color="black",
                   
                   point.padding = .1,
                   #family="Open Sans",
                   alpha=.2,
                   
                   size=3,
                   
                   #size=2,
                   min.segment.length = 0,
                   segment.alpha=.25,
                   max.overlaps = 500000000,
                   force = 5,
                   hjust= 0.5,
                   vjust=0.5,
                   box.padding = 0,
                   show.legend = FALSE)+
  
  
  geom_hline(yintercept = 5, alpha=.6, linetype='dotted')+
  
  
  scale_x_continuous(
    breaks = breaks_dynamic,  # Breaks at each integer position
    limits = limits_dynamic,    # Padding on both ends
    labels = pathway_list_factor$pathway_list # Custom labels
  )+
  
  labs(x="pathway term group", 
      y="-log10(p value)",
      title = str_wrap(names(gost_v1)[[de2]], width = 30),
      intersection_size="Intersection Size")+

  guides(fill=guide_legend(title="Pathway Term Group"),
         size=guide_legend(title="Intersection Size"))+
  
  continuous_scale(
    aesthetics = c("size", "point.size"),  # Scale both point size and label size
    palette = my_pal(c(1, 8))  # Custom palette for scaling
  ) +
  coord_cartesian(clip = "off")+
  theme_classic()





# 
# 
# 
# plotinloop<-gost_v1[[de2]]$plot1
# 
# 
# built <- ggplot_build(plotinloop)
# label_data <- built$data[[3]]  # contains x, y, label positions
# 
# 
# 
# 
# 
# 
# plotinloop+
# 
# geom_label(data = label_data,
#            aes(x = x, y = y, label = label),
#            color = "black", size = 3) 
# 
# 







#print(gost_v1[[de2]]$plot)

#get generation output metadata.
metadata_text <- paste0("src: ", 
                        rstudioapi::getSourceEditorContext()$path %>%
                          sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.) %>%
                          sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.), 
                        "\n", "time: ", 
                        lubridate::round_date(Sys.time(), "second"))

text_grob <- textGrob(metadata_text, x=.05, just="left", gp = gpar(fontsize = 7, col = "gray30"))

#add in p value annotation.
gost_v1[[de2]]$plot2 <- arrangeGrob(
  gost_v1[[de2]]$plot1, 
  text_grob, 
  ncol = 1, 
  heights = c(3, 0.3),  # Maintain spacing
  layout_matrix = rbind(c(1), c(2))
)


}














# grid.newpage()
# grid.draw(gost_v1[[7]]$plot1)
# 
# grid.newpage()
# grid.draw(gost_v1[[7]]$plot2)
# gost_v1[[7]]$result%>%as_tibble()%>%
#   view()


plot_output_base<-"gost-"
plot_output_loc<-"nf1g/gost/plots/"  #include trailing slash


# Loop through the outer list
for (ab in seq_along(names(gost_v1))) {
  
  x1<-names(gost_v1)[ab]
  
  xtitle<-x1%>%gsub(".rds", "", .)
  
  # cat("1: ", red(x1), "\n")
  # 
  # cat("2: ", blue(xtitle), "\n")
  # 

  # Generate the filename dynamically
  file_name <- fs::path_sanitize(paste0(
      plot_output_base,
      xtitle,
      ".pdf"
    ))
# 
#   cat("3: ", green(file_name), "\n")
  
  plot_object<-gost_v1[[ab]]$plot2
  
  grid.newpage()
  #grid.draw(plot_object)
  
  
  title_text <- paste0(

    gost_v1[[ab]]$plot1$labels$title %>%
      gsub("\n", "][", .) %>%
      gsub(".rds", "", .)%>%
      gsub(";", "", .) %>%
      gsub("\\,", "", .),

    "src: ",
    rstudioapi::getSourceEditorContext()$path %>%
      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/", "", .) %>%
      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/", "", .),
    " at ", round_date(Sys.time(), "second")
  )

  

  # showtext::showtext_opts(dpi = 600)
  ggsave(
    path = plot_output_loc,
    filename = file_name,
    plot = plot_object,
    title = title_text,
    limitsize = FALSE,
    height = 7,
    width = 6,
    scale = 1.2,
    dpi = 600
  )
  
 

    # Print status message
    cat(bold(green("Saved plot: ")), cyan(file_name, "\n"), blue("in "), plot_output_loc, "\n")

}






plot_pdf_paths <-fs::dir_info(plot_output_loc, 
                              pattern = paste0("^", plot_output_base, ".*\\.pdf$", recurse = FALSE)) %>%
  filter(modification_time>Sys.time()-lubridate::minutes(40))%>%
  filter(type == "file") %>%
  filter(!grepl("^comb", basename(path)))%>%
  arrange(modification_time)%>%
  tibble()%>%
  pull(path)

qpdf::pdf_combine(plot_pdf_paths, 
                  output = paste0(plot_output_loc,  "comb-gost-all.pdf"))






ggsave(
  path = plot_output_loc,
  filename = file_name,
  plot = grid.draw(gost_v1[[ab]]$plot2),
  title = title_text,
  limitsize = FALSE,
  height = 5,
  width = 10,
  scale = 1.2,
  dpi = 600
)

