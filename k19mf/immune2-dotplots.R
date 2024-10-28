source("libs.R")
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)
library(patchwork)


###########################################
all_sample_info<-readRDS("k19mf/ds/vm-00-sample_info.rds")

rpkms<-readRDS("k19mf/ds/immune2_rpkms.rds")

stats<-readRDS("k19mf/ds/immune2-acral-v-subq-stats.rds")

groups<-readRDS("k19mf/ds/immune2-genes_groups.rds")



immune_data<-rpkms%>%select(gene_name_ms, sample_id, rpkm)%>%
  left_join(all_sample_info%>%select(mouse_num, sample_id, tumor_type))%>%
  left_join(stats)










# 
# gene_names<-int2%>%select(gene_name_ms)%>%unique()%>%pull()

groups_list <- groups %>%
  group_by(group) %>%
  summarise(gene_name_ms = list(gene_name_ms)) %>%
  deframe()



group_names<-names(groups_list)
# group_names<-names(groups_list)[1]


for (gr in 1:length(group_names)){

group_name<-group_names[gr]


genes_in_group<-groups_list[[group_name]]



plot_list<-list()
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
  
  plot_list[[ge]]<-p
    
}else{

  plot_list[[ge]]<-p+
    geom_signif(comparisons = list(c("acral", "subq")),
                margin_top = 0.08,
                tip_length = 0,
                textsize = 6,
                annotations = c(int_sub$fdr_stars[1]))  # Specify your p-value here
    #y_position = 9)+

}
}



# plot_list
# 

plot_grid<-grid.arrange(grobs=plot_list, ncol=3,
                       
                      top = textGrob(group_name, gp = gpar(fontsize = 20, fontface = "bold")))




quotient <- length(plot_list) %/% 3  # Integer division, gives 4
remainder <- length(plot_list) %% 3  # Modulus, gives 1

if (remainder == 0){
  
height_rows=quotient  
  
  
}else{
  
height_rows=quotient+1  

}




ggsave(paste0("k19mf/plots/", "immune2-dotplots-", group_name, ".pdf"),
       plot=plot_grid,
       height=6*height_rows+1,
       width=10,
       scale = .75,
       dpi=600,
       limitsize = FALSE

)



}















plots<-fs::dir_info("k19mf/plots", full.names = T, pattern="^immune.*\\.pdf$", recurse = F)%>%tibble()%>%
  filter(modification_time>Sys.time()-lubridate::minutes(5))%>%
  arrange(path)%>%
  filter(type=="file")%>%
  pull(path)



qpdf::pdf_combine(plots, output = "k19mf/plots/combined/immune2-dotplots.pdf")

