source("libs.R")
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)



###########################################
sample_info<-readRDS("acral_paired/ds/v00-sample_info.rds")

rpkms<-readRDS("acral_paired/ds/v02-filtered_rpkms.rds")


gene_sigs<-readRDS("acral_paired/dexps/dexp-pairedsample_type-acral v. subq.rds")%>%
  select(gene_name_ms, PValue)



goi0<-rpkms%>%group_by(gene_name_ms)%>%
  slice(1)%>%
  select(1:4)%>%
  ungroup()%>%
  
  left_join(gene_sigs)

genes_list<-c("crkl", "nf1", "gab2", "tert", "mage", "gpr19", "mmp14", "rreb1")

goi1<-goi0%>%
  filter(grepl("^hox", gene_name_ms, ignore.case=TRUE)|
           grepl("^bcl2a1", gene_name_ms, ignore.case=TRUE)|
           grepl("^myb", gene_name_ms, ignore.case=TRUE)|
           grepl("^cep6", gene_name_ms, ignore.case=TRUE)|
           tolower(gene_name_ms) %in% genes_list)




df<-rpkms%>%
  
  semi_join(goi1, by="gene_name_ms")%>%
  left_join(sample_info%>%select(mouse_num, sample_type, sample_id))


goi<-goi1%>%pull(gene_name_ms)

# goi<-goi[1]


plot_list<-list()

for (g in 1:length(goi)){

  
  

  
gene<-goi[g]  


pval<-gene_sigs%>%filter(gene_name_ms==!!gene)%>%pull("PValue")
  

df1<-df%>%
  filter(tolower(gene_name_ms)==tolower(gene))

x_label <- unique(df1$gene_name_ms)[1]


plot_list[[g]]<-ggplot(df1, aes(x=sample_type, y=rpkm))+
  geom_line(aes(group = mouse_num), alpha=.3)+
  geom_point(aes(color=sample_type),size=4, alpha=.7)+
  geom_point(size=4, shape=21, alpha=.8)+
  geom_text(aes(label=ifelse(sample_type=="acral",
                             as.character(mouse_num),'')),hjust=1.5,vjust=0.5,
            size = 3)+

  
  # 
  # facet_grid("gene_name_ms")+
  scale_color_manual(values=c("subq"="#fffa99", "acral"= "#99fff6"))+
  labs(title=x_label,
       x="")+
  annotate(
    "text", 
    label = paste0("p = ", round(pval, 6)),
    x = 2, y = 15, size = 5, colour = "black"
  )+
 
  theme_classic()+
theme(legend.position = "none",
      aspect.ratio=2)

}


# Combine the plots into a grid
combined_plot <- grid.arrange(grobs = plot_list)









ggsave(paste0("acral_paired/plots/", "combined_individual_genes_dot-paired",".pdf"),
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       
       
       
       
       plot=combined_plot,
       
       
       
       
       limitsize = FALSE,
       
       
       height=40,
       width=40,
       scale = 1,
       dpi=600,
       
       
       
)

stop()
































genes_to_keep<-stats%>%
  filter(fdr_stars!="ns")%>%
  group_by(gene_name_ms)%>%slice(1)%>%ungroup()%>%
  pull(gene_name_ms)





immune_data<-rpkms%>%select(gene_name_ms, sample_id, rpkm)%>%
  left_join(all_sample_info%>%select(mouse_num, sample_id, tumor_type))%>%
  left_join(stats)%>%
  
  filter(gene_name_ms %in% genes_to_keep)



#subset groups - gene pairs by significance of genes. 
groups2<-groups%>%
  filter(gene_name_ms %in% genes_to_keep)%>%
  select(group, gene_name_ms)

#turn group - gene pairs into named list. 
groups_list <- groups2 %>%
  group_by(group) %>%
  summarise(gene_name_ms = list(gene_name_ms)) %>%
  deframe()


#names of groups that are included.
group_names<-names(groups_list)


plot_list<-list()
for (gr in 1:length(group_names)){
  
  group_name<-group_names[gr]
  
  
  genes_in_group<-groups_list[[group_name]]
  
  
  gene_plots <- list()
  
  for(ge in 1:length(genes_in_group)){
    
    
    gene_name<-genes_in_group[ge]
    
    
    int_sub<-immune_data%>%filter(gene_name_ms==gene_name)
    
    
    p<-ggplot(int_sub, aes(x=tumor_type, y=rpkm))+
      
      geom_point(aes(color=tumor_type),size=4, alpha=.7)+
      geom_point(size=4, shape=21, alpha=.8)+
      
      scale_color_manual(values=c("subq"="#fffa99", "acral"= "#99fff6"))+
      
      
      geom_errorbar(aes(y=mean, ymin=mean-se, ymax=mean+se),
                    width=.07, linewidth=.3)+
      stat_summary(aes(y = rpkm, ymax = after_stat(y), ymin = after_stat(y)),
                   fun = mean, geom = "errorbar", color = "black", width = 0.2, linewidth=.3)+
      
      
      
      scale_y_continuous(
        limits = c(NA, max(int_sub$rpkm)*1.1),                       # Allow full range of y values
        breaks = seq(0, max(int_sub$rpkm)*1.05, by = round(max(int_sub$rpkm)+.5)/5)              # Specify where tick marks appear
      )+
      
      
      labs(title = int_sub$gene_name_ms[1],
           x="Tumor Type",
           y="RPKM")+
      
      theme_classic()+
      theme(legend.position = "none",
            plot.title = element_text(size = 18),
            #axis.line.y = element_blank(),
            plot.margin = margin(c(5,5,5,5))
      )
    
    
    
    
    if(int_sub$fdr_stars[1]=="ns"){
      
      gene_plots[[ge]]<-p
      
    }else{
      
      gene_plots[[ge]]<-p+
        geom_signif(comparisons = list(c("acral", "subq")),
                    margin_top = 0.08,
                    tip_length = 0,
                    textsize = 6,
                    annotations = c(int_sub$fdr_stars[1]))  # Specify your p-value here
      #y_position = 9)+
      
    }
  }
  
  plot_list[[group_name]] <- gene_plots
  
  
}




library(cowplot)




ncol<-max(sapply(plot_list, length))+1





# Create a list to hold the plot grids for each gene group
plot_grid_list <- lapply(plot_list, function(group_plots) {
  
  
  # Pad the group with NULLs if it has fewer than `ncol` plots
  padded_group_plots <- c(group_plots, rep(list(NULL), ncol - length(group_plots) %% ncol))
  
  # Arrange each gene group in its own row, maintaining the same number of columns
  plot_grid <- plot_grid(plotlist = padded_group_plots, ncol = ncol, align = "hv")
  
  return(plot_grid)
  
  
})


plot_grid_list <- lapply(names(plot_list), function(group_name) {
  group_plots <- plot_list[[group_name]]
  
  # Pad the group with NULLs if it has fewer than `ncol` plots
  padded_group_plots <- c(group_plots, rep(list(NULL), ncol - length(group_plots) %% ncol))
  
  # Arrange the plots in a grid
  plot_grid <- plot_grid(plotlist = padded_group_plots, ncol = ncol, align = "hv")
  
  # Create a title for the group
  title_plot <- ggdraw() + 
    draw_label(group_name, fontface = 'bold', size = 14, x = 0.5, y = 0.5, hjust = 0.5)
  
  # Combine the title and the plot grid
  combined_plot <- plot_grid(title_plot, plot_grid, ncol = 1, rel_heights = c(0.1, 1))
  return(combined_plot)
})



plot_grid_list <- lapply(names(plot_list), function(group_name) {
  group_plots <- plot_list[[group_name]]
  
  # Pad the group with NULLs if it has fewer than `ncol` plots
  padded_group_plots <- c(group_plots, rep(list(NULL), ncol - length(group_plots) %% ncol))
  
  # Arrange the plots in a grid
  plot_grid <- plot_grid(plotlist = padded_group_plots, ncol = ncol, align = "hv")
  
  # Create a title for the group
  title_plot <- ggdraw() + 
    draw_label(group_name, fontface = 'bold', size = 10, x = .5, y = 0.5, hjust = .5, vjust = 0.5)
  
  # Combine the title and the plot grid side by side
  combined_plot <- plot_grid(title_plot, plot_grid, ncol = 2, rel_widths = c(0.2, 1))
  return(combined_plot)
})

# Set equal row heights for each gene group to avoid squishing
final_plot <- plot_grid(plotlist = plot_grid_list, ncol = 1, align = "v", rel_heights = rep(1, length(plot_grid_list)))




# Display the final plot
print(final_plot)



ggsave(paste0("k19mf/plots/", "immune2-dotplots-sigs-",  ".pdf"),
       plot=final_plot,
       height=35,
       width=35,
       scale = .75,
       dpi=600,
       limitsize = FALSE
       
)



# 
# 
# quotient <- length(plot_list) %/% 3  # Integer division, gives 4
# remainder <- length(plot_list) %% 3  # Modulus, gives 1
# 
# if (remainder == 0){
#   
# height_rows=quotient  
#   
#   
# }else{
#   
# height_rows=quotient+1  
# 
# }
# 
# 
# 
# 

# 
# 
# 
# 
# 
# 
# 












plots<-fs::dir_info("k19mf/plots", full.names = T, pattern="^immune.*\\.pdf$", recurse = F)%>%tibble()%>%
  filter(modification_time>Sys.time()-lubridate::minutes(5))%>%
  arrange(path)%>%
  filter(type=="file")%>%
  pull(path)



qpdf::pdf_combine(plots, output = "k19mf/plots/combined/immune2-dotplots-sig.pdf")








