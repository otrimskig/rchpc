source("libs.R")
library(tidyverse)
library(ggplot2)
library(EnhancedVolcano)

df<-readRDS("nf1g/tumor_types/aod-gene-analysis.rds")




ggplot(df, aes(x=-log(hazard_ratio), y=-log(p_value)))+
  geom_point()



# Example thresholds
pval_threshold <- 0.01
hr_threshold <- log(2) 

df$significance <- factor(
  ifelse(abs(log(df$hazard_ratio)) > hr_threshold & df$p_value < pval_threshold,
         "Significant", "Not Significant")
)

# Expand axis limits
x_limits <- range(df%>%filter(hazard_ratio!=0)%>%pull(hazard_ratio)%>%-log(.), na.rm = TRUE) * 2000
y_limits <- range(-log10(df$p_value), na.rm = TRUE) * 1.1

p<-ggplot(df, aes(x = -log(hazard_ratio), y = -log10(p_value), color = significance, label = gene)) +
  geom_point(alpha = 0.7)+
  geom_text_repel(data = subset(df, significance == "Significant"), max.overlaps = 5, color="black")+
  
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "purple"))+
  
  
  geom_vline(xintercept = c(-hr_threshold, hr_threshold), alpha=.5)+
  geom_hline(yintercept = -log10(pval_threshold)) +
  
  
  labs(x = "-log(Hazard Ratio)", y = "-log10(P-value)", color = "Significance")+
# Add 10% padding on both sides
 
 # Expand axis limits
  
  theme_minimal()



















ggsave("nf1g/tumor_types/vol-haz.pdf",
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       plot=p,
       limitsize = FALSE,
       
       
       height=10,
       width=10,
       scale = 1,
       dpi=600,
       
       
       
)
















stop()



  

dfl<-df%>%
  mutate(hazard_ratio=-log(hazard_ratio))


EnhancedVolcano(dfl, 
                   lab=dfl$gene,
                   x="hazard_ratio",
                   y="p_value",
                   axisLabSize = 12,
                   labSize = 5.0,
                   drawConnectors = TRUE,
                   boxedLabels = FALSE,
                   colConnectors = "grey50",
                   arrowheads = FALSE,
                   min.segment.length=.3,
                   #title=plot_title,
                   subtitle = "Of major category: Gliomas",
                   maxoverlaps = 100,
                   #typeConnectors="open",
                   
                   #ylim=c(0, y_axis_max),
                   #xlim=c(axes$minfc,axes$maxfc),
                   
                   endsConnectors="first",
                   legendPosition = 'right',
                   legendLabSize = 10,
                   legendIconSize = 4.0,
                   col = c('grey', 'light blue', 'light green', 'purple')
                   
                   
                   #selectLab = c("Pten", "Cdkn2a", "Nf1", "Atrx")
                   
)

