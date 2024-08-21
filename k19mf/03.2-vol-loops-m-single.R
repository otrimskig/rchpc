
source("libs.R")
library(tidyverse)
library(EnhancedVolcano)


###########################################

#get vector of all .rds files containing dexp analyses.
#all_dexps<-list.files("m-dexps", full.names = T, pattern="^dexp.*\\.rds$")
all_dexps<-"k19mf/dexps/dexp-tumor_type-subq v. foot padexcluded2s.rds"

all_dexps<-"k19mf/dexps/dexp-tumor_type-subq v. foot pad.rds"

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



#now load the actual dataset from the predetermined path.
de_df<-readRDS(rds_file_path)

vdf<-de_df


#axes<-readRDS("ds/de_axis_lims.rds")


minpv<-vdf%>%
  summarise(pv=min(PValue))%>%
  pull()




y_axis_max<-max(-log10(minpv),  12)


p<-EnhancedVolcano(vdf, 
                lab=vdf$gene_name_ms,
                
                
                x="logFC",
                y="PValue",
                axisLabSize = 20,
                labSize = 8.0,
                
                #drawConnectors = TRUE,
                
                
                boxedLabels = FALSE,
                colConnectors = "grey50",
                arrowheads = FALSE,
                
                min.segment.length=1,
               
                
                title=plot_title,
                #subtitle = "Of major category: Gliomas",
                maxoverlaps = 3,
                #typeConnectors="open",
                
                #ylim=c(0, y_axis_max),
                #xlim=c(axes$minfc,axes$maxfc),
                
                endsConnectors="first",
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 15.0,
                col = c('grey', 'light blue', 'light green', 'purple')
                
                
                #selectLab = c("Pten", "Cdkn2a", "Nf1", "Atrx")
                
                )


#pdf_height<-y_axis_max/15+2
pdf_height=8

ggsave(paste0("k19mf/plots/",fs::path_sanitize(paste0("vol-glioma-atrx", file_base, "-rpkms.pdf"))),
       plot=p,

       scale = 1,
       dpi=600,
       width = 20,
       height =20,
       unit="in",
       limitsize = FALSE

)


}




stop("continue if you need to combine all pdfs.")

plots<-list.files("k19mf/plots", full.names = T, pattern="^vol.*\\.pdf$")

plots<-list.files("k19mf/plots", full.names = T, pattern="\\.pdf$")



qpdf::pdf_combine(plots, output = "k19mf/plots/combined/all.pdf")
