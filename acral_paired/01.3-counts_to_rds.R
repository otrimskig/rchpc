library(tidyverse)
library(dtplyr)
library(data.table)

#set directory of output file (relative to rchpc proj)
output_path<-"acral_paired/ds/"

exp_path<-"../exp_data/25622R"

counts_files<-c(list.files(path="../exp_data/25622R/featurecounts",  pattern="FeatureCounts.txt$", full.names = TRUE))




file_path<-counts_files[1]

counts<-fread(file_path)

# counts[,3]


all_counts<-counts[,1:2]


for (i in 1:length(counts_files)){
  
file_path<-counts_files[i]

counts<-fread(file_path)

all_counts<-cbind(all_counts, counts[,3])

}


saveRDS(all_counts, paste0(output_path, "v00-all_counts.rds"))

