if (!exists("n.cores")) {
  # Run your code here
"initilizing cores..."
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

"parallel cores initialized."
}


library(tidyverse)
# axes<-readRDS("ds/de_axis_lims.rds")
# axes$pv

# -log10(axes$pv)
# 
# -log10(axes$pv)


###########################################

#get vector of all .rds files containing dexp analyses.
all_dexps<-list.files("m-dexps", full.names = T, pattern="^dexp.*\\.rds$")


library(foreach)

foreach(i=1:length(all_dexps)) %dopar% {

library(tidyverse)
library(EnhancedVolcano)

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


axes<-readRDS("ds/de_axis_lims.rds")


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
                subtitle = "",
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

ggsave(paste0("plots/",fs::path_sanitize(paste0("vol-", file_base, "-rpkms.pdf"))),
       plot=p,

       scale = 1.2,
       dpi=600,
       width = 10,
       height =pdf_height,
       unit="in"

)


}



plots<-list.files("plots", full.names = T, pattern="^vol.*\\.pdf$")



qpdf::pdf_combine(plots, output = "plots/combined/vol-scale_same.pdf")
