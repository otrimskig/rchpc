source("libs.R")
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gridExtra)

int<-readRDS("k19mf/ds/int-pvals.rds")


int_sum<-int%>%
  group_by(gene_name_ms, tumor_type)%>%
  summarise(mean=mean(rpkm),
            sd=sd(rpkm),
            n=n(),
            se=sd/sqrt(n))


int2<-int%>%left_join(int_sum)%>%
  mutate(stars=gtools::stars.pval(p_adj_fdr))%>%
  mutate(stars=trimws(gsub("\\.", "", stars)))%>%
  mutate(stars=gsub("^$", "ns", stars))



gene_names<-int2%>%select(gene_name_ms)%>%unique()%>%pull()



plot_list<-list()



for (g in 1:length(gene_names)){
gene_name<-gene_names[g]

int_sub<-int2%>%filter(gene_name_ms==gene_name)


if(int_sub$stars[1]=="ns"){
  
  plot_list[[g]]<-ggplot(int_sub, aes(x=tumor_type, y=rpkm))+
    
    
    geom_point(aes(color=tumor_type),size=4, alpha=.7)+
    geom_point(size=4, shape=21, alpha=.8)+
    
    scale_color_manual(values=c("subq"="#fffa99", "acral"= "#99fff6"))+
    
    
    
    
    
    geom_errorbar(aes(y=mean, ymin=mean-se, ymax=mean+se),
                  width=.07, size=.3)+
    stat_summary(aes(y = rpkm, ymax = after_stat(y), ymin = after_stat(y)),
                 fun = mean, geom = "errorbar", color = "black", width = 0.2, size=.3)+
    
    # stat_pvalue_manual(
    #   stat.test, x = "tumor_type", y.position = 33,
    #   label = "p.signif",
    #   position = position_dodge(0.8)
    # )+
    
    # geom_signif(comparisons = list(c("acral", "subq")), 
    #             margin_top = 0.08,
    #             tip_length = 0,
    #             textsize = 6,
    #             annotations = c(int_sub$stars[1]))+  # Specify your p-value here
    # #y_position = 9)+
    
    scale_y_continuous(
      limits = c(NA, max(int_sub$rpkm)*1.1),                       # Allow full range of y values
      breaks = seq(0, max(int_sub$rpkm)*1.05, by = round(max(int_sub$rpkm)/5))              # Specify where tick marks appear
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
    
    
    # annotation_custom(
    #   grid::linesGrob(gp = grid::gpar(lwd = 3, col = "black")), 
    #   xmin = -Inf, xmax = -Inf, ymin = -4, round(max(int_sub$rpkm)*1.05)  # Control the height of the y-axis line
    # )
  
  
}else{


plot_list[[g]]<-ggplot(int_sub, aes(x=tumor_type, y=rpkm))+
 

  geom_point(aes(color=tumor_type),size=4, alpha=.7)+
  geom_point(size=4, shape=21, alpha=.8)+
  
  scale_color_manual(values=c("subq"="#fffa99", "acral"= "#99fff6"))+
  
  
  
  
  
  geom_errorbar(aes(y=mean, ymin=mean-se, ymax=mean+se),
                width=.07, size=.3)+
  stat_summary(aes(y = rpkm, ymax = after_stat(y), ymin = after_stat(y)),
               fun = mean, geom = "errorbar", color = "black", width = 0.2, size=.3)+
  
  # stat_pvalue_manual(
  #   stat.test, x = "tumor_type", y.position = 33,
  #   label = "p.signif",
  #   position = position_dodge(0.8)
  # )+
  
  geom_signif(comparisons = list(c("acral", "subq")), 
              margin_top = 0.08,
              tip_length = 0,
              textsize = 6,
              annotations = c(int_sub$stars[1]))+  # Specify your p-value here
              #y_position = 9)+
 
  scale_y_continuous(
    limits = c(NA, max(int_sub$rpkm)*1.1),                       # Allow full range of y values
    breaks = seq(0, max(int_sub$rpkm)*1.05, by = round(max(int_sub$rpkm)/5))              # Specify where tick marks appear
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
  
  
  # annotation_custom(
  #   grid::linesGrob(gp = grid::gpar(lwd = 3, col = "black")), 
  #   xmin = -Inf, xmax = -Inf, ymin = -4, ymax = round(max(int_sub$rpkm)*1.05)  # Control the height of the y-axis line
  # )

}

}


plot_list[[1]]



plot_grid<-grid.arrange(grobs=plot_list, ncol=3)

#plot_grid



ggsave(paste0("k19mf/plots/", "hm-integrins-dotplots.pdf"),
       plot=plot_grid,
       height=18,
       width=4,
       scale = 2.5,
       dpi=600,
       limitsize = FALSE
       
)




df<-`int-pvals`

df2<-df%>%group_by(gene_name_ms)%>%slice(1)%>%
  select(gene_name_ms, starts_with("p_"))

write_csv(df, "k19mf/ds/int_pvals.csv")
write_csv(df2, "k19mf/ds/int_pvals2.csv")
