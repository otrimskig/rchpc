source("libs.R")

library(tidyverse)
library(readxl)
library(janitor)
library(purrr)

library(foreach)

if (!exists("n.cores")) {
  
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




files<-fs::dir_ls("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/rchpc/nec_files/nec")%>%
  as_tibble()


files2<-files%>%
  rename(path=value)%>%
  mutate(fn=basename(path))%>%
  mutate(ft=if_else(grepl("xlsx$", fn), "xlsx", 
                    if_else(grepl("xls$", fn), "xls", NA)))%>%
  filter(ft=="xls"|ft=="xlsx")




all_sheets<-list()

for (i in 1:length(files2$path)){


filename<-basename(files2$path[i])

sheet_names<-excel_sheets(files2$path[i])


for (s in 1:length(sheet_names)){

all_sheets[[filename]][[s]]<-sheet_names[s]
  
}


}



fl<-tibble(
  file = rep(names(all_sheets), times = sapply(all_sheets, length)),
  sheet = unlist(all_sheets)
)%>%
  
  mutate(file_sheet=paste(file, sheet))%>%
  mutate(path=paste0("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/rchpc/nec_files/nec/", file))








all_files_list<-list()






x<-foreach(i=1:nrow(fl)) %dopar% {
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(tidyr)




read_excel(fl$path[i], sheet = fl$sheet[i]) %>%
    clean_names() %>%
    #mutate(across(everything(), ~ ifelse(. == "`", NA, .))) %>%
    select(where(~ !all(is.na(.)))) %>%
    t() %>%
    as_tibble() %>%
    row_to_names(1) %>%
    clean_names() %>%
    mutate(filename = fl$file[i])%>%
    mutate(sheet_name2 = as.character(fl$sheet[i]))


}





saveRDS(x, "nec_files/outside_files.rds")

