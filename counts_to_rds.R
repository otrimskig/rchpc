library(tidyverse)
library(dtplyr)
library(data.table)



counts_files<-list.files(path="../exp_data/23908R_merged/featurecounts/",  pattern="FeatureCounts.txt$", full.names = TRUE)
file_path<-counts_files[1]

counts<-fread(file_path)

counts[,3]


all_counts<-counts[,1:2]


for (i in 1:length(counts_files)){
  
file_path<-counts_files[i]

counts<-fread(file_path)

all_counts<-cbind(all_counts, counts[,3])

}

