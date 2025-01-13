##########################
source("libs.R")
library(tidyverse)

sample_info<-readRDS("acral_paired/ds/v00-sample_info.rds")

sample_name_map<-sample_info%>%
  select(sample_id, mouse_num, sample_type)%>%
  mutate(mouse_num_type_clean=paste0("x", mouse_num, "_", sample_type))


saveRDS(sample_name_map, "acral_paired/ds/name_map.rds")


#add dexp set or other df needed to be renamed.
df<-readRDS(NULL)


df2<-df%>%
  rename_with(~ {
    # For each column name, apply string replacements based on the mapping
    updated_names <- .
    for (i in seq_along(sample_name_map$sample_id)) {
      updated_names <- str_replace_all(
        updated_names,
        fixed(sample_name_map$sample_id[i]),
        sample_name_map$mouse_num_type_clean[i]
      )
    }
    updated_names
  })



