library(tidyverse)


neds<-c(list.files("m-dexps", pattern="NED", full.names = T), 
        list.files("m-dexps", pattern="No evidence of disease", full.names = T))


neds_filtered<-tibble(file_path1=neds)%>%
  mutate(basename1=basename(file_path1))%>%
  
  mutate(comp=sub("dexp-","",basename1))%>%
  mutate(category= sub("-.*", "", comp))%>%
  mutate(category= sub("-.*", "", comp))%>%
  
  mutate(groups=sub(".*-", "", comp))%>%

  mutate(group1=sub(" v\\..*", "", groups))%>%
  mutate(group2=sub(".* v\\. ", "", groups))%>%
  mutate(group2=sub("\\.rds", "", group2))%>%
  
  
  filter(group1!="NED")%>%
  filter(group1!="No evidence of disease")%>%
  
  mutate(file_path2=paste0("m-dexps/", "dexp-", category, "-", group2, " v. ", group1, ".rds"))%>%
  
  mutate(file_path3=paste0("m-dexps/switched/", basename1))%>%
  
  relocate(file_path1, file_path2, file_path3)




for(i in 1:length(neds_filtered$file_path1)){

path1<-neds_filtered$file_path1[i]
path2<-neds_filtered$file_path2[i]
path3<-neds_filtered$file_path3[i]

readRDS(path1)%>%
  mutate(logFC=logFC*-1)%>%
saveRDS(path2)

fs::file_move(path1, path3)
}
