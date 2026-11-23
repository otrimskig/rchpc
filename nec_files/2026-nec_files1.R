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




files<-fs::dir_ls("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/necfiles0/files0")%>%
  as_tibble()


files2<-files%>%
  rename(path=value)%>%
  mutate(fn=basename(path))%>%
  mutate(ft=if_else(grepl("xlsx$", fn), "xlsx", 
                    if_else(grepl("xls$", fn), "xls", NA)))



all_files_list<-list()



x<-foreach(i=1:length(files2$path)) %dopar% {
  source("libs.R")
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(tidyr)

# for(i in 1:length(files2$path)){

read_excel(files2$path[i])%>%
  clean_names()%>%
  mutate(across(everything(), ~ ifelse(. == "`", NA, .)))%>%
  select(where(~ !all(is.na(.))))%>%
  t()%>%
  as_tibble()%>%
  row_to_names(1)%>%
  clean_names()%>%
  
  mutate(filename=files2$fn[i])


}


# Stop the parallel cluster when you're done
if (exists("my.cluster")) {
  parallel::stopCluster(my.cluster)
  rm(my.cluster)  # Optional: remove the cluster object from environment
  message("Parallel cluster stopped. Returning to single-core execution.")
}







x1 <- keep(x, ~ nrow(.) >= 1)


x2 <- bind_rows(x1)


x2.1<-x2%>%select(where(~ !all(is.na(.))))%>%
  
  relocate(filename)%>%
  unite(
    disposition_notes,
    starts_with("dispos"),
    sep = "; ",
    na.rm = TRUE,
    remove = TRUE
  )%>%
  mutate(disposition_notes=if_else(disposition_notes=="", NA, disposition_notes))%>%

unite(
  additional_comments,
  starts_with("additional_com"),
  sep = "; ",
  na.rm = TRUE,
  remove = TRUE
)%>%
  mutate(additional_comments=if_else(additional_comments=="", NA, additional_comments))%>%
  
  
  mutate(across(everything(), ~ {
    # Apply logic to each value in the column
    sapply(., function(x) {
      # Check if the value is a 5-digit number between 40000 and 46000
      if (grepl("^\\d{5}$", as.character(x)) && as.numeric(x) >= 40000 && as.numeric(x) <= 46000) {
        # Convert the value into an Excel date and then back to character
        return(as.character(as.Date(as.numeric(x), origin = "1899-12-30")))
      }
      # Return the value as-is if it doesn't meet the condition
      return(as.character(x))
    })
  }))%>%
  
  
  
  mutate(genotype=if_else(is.na(genotype)&grepl("^Nes", na),  na, genotype))%>%
  mutate(na=if_else(genotype==na, NA, genotype))





saveRDS(x2.1, "nec_files/26nec0.rds")
