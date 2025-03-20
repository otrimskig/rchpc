source("libs.R")

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(gprofiler2)
library(ggpp)
library(crayon)
library(grid)
library(gridExtra)

#get all dexp files.
dexp_files <- list.files(
  path = "nf1g/dexps/atrx_comps",
  pattern = "^dexp-.*\\.rds$",
  full.names = TRUE
)


#dexp_files<-dexp_files[1]

gost_v1<-list()
for (de1 in 1:length(dexp_files)){
  
  #de1 <- 1  # Ensure i starts at 1 for debugging
  
  # Load and filter data
  dexp <- readRDS(dexp_files[de1]) %>%
    filter(FDR < 0.1) %>%
    filter(abs(logFC) > 1) %>%
    arrange(logFC)
  
  dexp_name <- basename(dexp_files[de1])
  
  cat(blue("working on "), red(dexp_name), " (", cyan(de1), " of ", length(dexp_files), ")", "\n")
  
  
  
  # Run gProfiler analysis
  gost_v1[[dexp_name]]<- gost(query = dexp$gene_id_ms,
            organism = "mmusculus",
            numeric_ns = "ENTREZGENE_ACC",
            ordered_query = TRUE)
  
  # # Extract results
  # results <- gdexp$result  # Extract enrichment results


}


saveRDS(gost_v1, "nf1g/atrx_comps/ds/gost_v1.rds")
