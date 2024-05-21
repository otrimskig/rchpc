library(tidyverse)
library(data.table)
library(dtplyr)
library(janitor)






# stat_file<-"../exp_data/23908R_merged/featurecounts/23908X01_merged.FeatureCounts_stats.txt"
# 
# df1<-fread(stat_file)
# 
# 
# sample_id<-basename(stat_file)%>%
#   substr(1,8)
# 
# t(df1)%>%as_tibble()%>%
#   row_to_names(1)%>%
#   janitor::clean_names()%>%
#  
#   mutate(across(everything(), trimws))%>%
#   mutate(across(everything(), as.numeric))%>%
# 
#   mutate(sample_id=sample_id)%>%
#   relocate(sample_id)



files<-list.files(path = "../exp_data/23908R_merged/featurecounts/", pattern = "stats.txt$", full.names = TRUE)
all_stats<-tibble()

for (i in 1:length(files)){

stat_file<-files[i]

df1<-fread(stat_file)


sample_id<-basename(stat_file)%>%
  substr(1,8)

bam_stat<-t(df1)%>%as_tibble()%>%
  row_to_names(1)%>%
  janitor::clean_names()%>%
  
  mutate(across(everything(), trimws))%>%
  mutate(across(everything(), as.numeric))%>%
  
  mutate(sample_id=sample_id)%>%
  relocate(sample_id)

all_stats<-rbind(all_stats, bam_stat)


}




combined_perc<-all_stats%>%
  mutate(total_reads = rowSums(select(., assigned:unassigned_ambiguity)))%>%
  mutate(across(matches("^assigned|^unassigned"), 
                ~ . / total_reads*100,
                .names = "{.col}_percent"))%>%
  select(sort(names(.)))%>%
  
  mutate(total_reads_m = round(total_reads/1E6, digits=2))%>%
  mutate(assigned_percent = round(assigned_percent, 2))%>%
  
  mutate(assigned_reads_m = round(assigned/1E6, digits=2))%>%
  relocate(sample_id, total_reads_m, assigned_percent, assigned_reads_m)


#
write_csv(combined_perc, file="../exp_data/23908R_merged/featurecounts/all_sample_stats.csv", na="")

  