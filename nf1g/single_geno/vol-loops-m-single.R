

library(tidyverse)
library(EnhancedVolcano)


###########################################

#get vector of all .rds files containing dexp analyses.
all_dexps<-list.files("nf1g/single_geno/dexps", full.names = T, pattern="^dexp.*\\.rds$")



for(i in 1:length(all_dexps)){



rds_file_path<-all_dexps[i]

#determine plot title from basic naming pattern of rds file.
plot_title<-basename(rds_file_path)%>%
  sub(".rds$", "", .)%>%
  sub("^dexp-", "", .)%>%
  sub(".*-", "", .)%>%
  sub(" v. ", " vs. ", .)


#get file base since we will need this later to actually name the file.
file_base<-basename(rds_file_path)%>%
  sub(".rds$", "", .)%>%
  sub("^dexp-", "", .)

genotype_analyzed<-sub("\\-.*", "", file_base)


#now load the actual dataset from the predetermined path.
de_df<-readRDS(rds_file_path)

vdf<-de_df


axes<-readRDS("nf1g/ds/de_axis_lims.rds")


minpv<-vdf%>%
  summarise(pv=min(PValue))%>%
  pull()




y_axis_max<-max(-log10(minpv),  12)


p<-EnhancedVolcano(vdf, 
                lab=vdf$gene_name_ms,
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
                
                
                
                subtitle = paste0("Of genotype: ", genotype_analyzed),
                
                
                
                maxoverlaps = 10,
                #typeConnectors="open",
                
                ylim=c(0, y_axis_max),
                xlim=c(axes$minfc,axes$maxfc),
                
                endsConnectors="first",
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 4.0,
                col = c('grey', 'light blue', 'light green', 'purple')
                
                
                #selectLab = c("Pten", "Cdkn2a", "Nf1", "Atrx")
                
                )


pdf_height<-y_axis_max/15+2
pdf_height=8

ggsave(paste0("nf1g/single_geno/plots/",fs::path_sanitize(paste0("vol-", file_base, "-rpkms.pdf"))),
       plot=p,
       
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/rchpc/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       

       scale = 1.2,
       dpi=600,
       width = 15,
       height =pdf_height,
       unit="in"

)


}

stop("stop here")

plots<-list.files("nf1g/single_geno/plots", full.names = T, pattern="^vol.*\\.pdf$")



qpdf::pdf_combine(plots, output = "nf1g/single_geno/plots/combined/vol-scale_same.pdf")
