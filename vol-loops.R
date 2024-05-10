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
library(EnhancedVolcano)
library(foreach)

###########################################
#set threshold values for deexps datasets.

# logfc_threshold<-2
# fdr_threshold<-.001

#get vector of all .rds files containing dexp analyses.
all_dexps<-list.files("dexps", full.names = T, pattern="^dexp.*\\.rds$")

#get sample metadata.
all_sample_info<-readRDS("ds/v07-per_sample_info.rds")

########################################
#volcano plot for same data.


foreach(i=1:length(all_dexps)) %dopar% {

library(tidyverse)
library(EnhancedVolcano)

rds_file_path<-all_dexps[i]

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
                
                endsConnectors="first",
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 4.0,
                col = c('grey', 'light blue', 'light green', 'purple')
                
                
                #selectLab = c("Pten", "Cdkn2a", "Nf1", "Atrx")
                
                )


ggsave(paste0("plots/",fs::path_sanitize(paste0("vol-", file_base, "-rpkms.pdf"))),
       plot=p,

       scale = 1.2,
       dpi=600,
       width = 10,
       height = 8,
       unit="in"

)


}



plots<-list.files("plots", full.names = T, pattern="^vol.*\\.pdf$")



qpdf::pdf_combine(plots, output = "plots/combined/cvol-1.pdf")
