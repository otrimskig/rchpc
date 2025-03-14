source("libs.R")

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


library(foreach)
library(tidyverse)
#library(gprofiler2)



#genes<-readRDS("ds/gene_stats.rds")


dexp_files<-list.files(path = "m-dexps", pattern= ".rds$", full.names = TRUE)





#sub("\\.rds$", ".pdf", dexp_files[1])
foreach(i=1:length(dexp_files)) %dopar% {
source("libs.R")
#parallel libraries.  
library(tidyverse)
library(gprofiler2)
# i<-1


html_out_path<-paste0("plots/gprof/", "gprof-", basename(sub("\\.rds$", ".html", dexp_files[i])))

d<-readRDS(dexp_files[i])%>%
  filter(FDR<0.01)%>%
  filter(abs(logFC)>2)%>%
  arrange(logFC)


d_name<-basename(dexp_files[i])


g<- gost(query = d$gene_id_ms,
               organism = "mmusculus",
               numeric_ns = "ENTREZGENE_ACC",
         ordered_query = TRUE)

p<-gostplot(g, capped = TRUE, interactive = T)




#save as self-contained html file. 
htmlwidgets::saveWidget(p, html_out_path)



}





# publish_gosttable(g)
