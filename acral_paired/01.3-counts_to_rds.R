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



#now clean up sample names. may need to adjust if format of names changes with other exp.


bam_sample_names<-colnames(all_counts)%>%as_tibble()%>%
  rename(name1=value)%>%
  mutate(name2=substr(name1, 1,9))

name_map<-setNames(bam_sample_names$name2, bam_sample_names$name1)

all_counts2<-all_counts%>%
  rename_with(~ name_map[.x], .cols = all_of(bam_sample_names$name1))



saveRDS(all_counts2, paste0(output_path, "v00-all_counts.rds"))


