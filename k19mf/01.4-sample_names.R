library(tidyverse)
library(dtplyr)
library(data.table)

#at this stage, check that the sample names (column titles on the counts matrix) are what you want them to be)

exp_dir<-"k19mf/"

all_counts<-readRDS("k19mf/ds/mv00-all_counts.rds")%>%head()

# tibble(a=colnames(all_counts)[3:ncol(all_counts)], b=NA)%>%
#   arrange(a)%>%
#   write_csv("k19mf/ds/a.csv")


b<-read_csv("k19mf/ds/b.csv")


rename_columns <- function(df, rename_table) {
  # Create a named vector for renaming
  rename_vector <- setNames(rename_table$b, rename_table$a)
  
  # Rename the columns in the data frame
  colnames(df) <- ifelse(colnames(df) %in% names(rename_vector),
                         rename_vector[colnames(df)],
                         colnames(df))
  
  return(df)
}

all_counts<-readRDS("k19mf/ds/mv00-all_counts.rds")

all_counts2 <- rename_columns(all_counts, b)%>%
  janitor::clean_names()


saveRDS(all_counts2, "k19mf/ds/vm-00-all_counts_sample_names_fixed.rds")

