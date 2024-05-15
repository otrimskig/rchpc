library(tidyverse)

#get vector of all .rds files containing dexp analyses.
all_dexps<-list.files("dexps", full.names = T, pattern="^dexp.*\\.rds$")

library(tidyverse)
library(EnhancedVolcano)


p<-30
  
  rds_file_path<-all_dexps[p]
  
  #determine plot title from basic naming pattern of rds file.
  plot_title<-basename(rds_file_path)%>%
    sub(".rds$", "", .)%>%
    sub("^dexp-", "", .)%>%
    sub("-", ": ", .)%>%
    sub("v", " vs. ", .)
  
  
  #get file base since we will need this later to actually name the file.
  file_base<-basename(rds_file_path)%>%
    sub(".rds$", "", .)%>%
    sub("^dexp-", "", .)
  
  
  
  #now load the actual dataset from the predetermined path.
  de_df<-readRDS(rds_file_path)
  
  vdf<-de_df
  
  
  axes<-readRDS("ds/de_axis_lims.rds")
  
  
  minpv<-vdf%>%
    summarise(pv=min(PValue))%>%
    pull()
  
  
  
  
  
y_axis_max<-round(max(-log10(minpv),  12)+1)
  

  
  p<-EnhancedVolcano(vdf, 
                     lab=vdf$gene_name,
                     x="logFC",
                     y="PValue",
                     axisLabSize = 12,
                     labSize = 5.0,
                     drawConnectors = TRUE,
                     boxedLabels = FALSE,
                     colConnectors = "grey50",
                     arrowheads = FALSE,
                     min.segment.length=.3,
                     title=plot_title,
                     subtitle = "",
                     maxoverlaps = 10,
                     #typeConnectors="open",
                     
                     ylim=c(-10, y_axis_max),
                     xlim=c(axes$minfc,axes$maxfc),
                     
                     endsConnectors="first",
                     legendPosition = 'right',
                     legendLabSize = 10,
                     legendIconSize = 4.0,
                     col = c('grey', 'light blue', 'light green', 'purple')
                     
                     
                     #selectLab = c("Pten", "Cdkn2a", "Nf1", "Atrx")
                     
  )
  
  
  
  
pdf_height<-(y_axis_max/20)+3
  
  ggsave("plots/test_plot7.pdf",
         plot=p,
         
         scale = 1.2,
         dpi=600,
         width = 10,
         height =pdf_height,
         unit="in"
         
  )
  
  