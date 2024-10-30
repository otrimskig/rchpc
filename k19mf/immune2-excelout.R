source("libs.R")
library(tidyverse)
library(dtplyr)
library(openxlsx)

stats<-readRDS("k19mf/ds/immune-acral-v-subq-stats.rds")


rpkms<-readRDS("k19mf/ds/immune2_rpkms.rds")%>%
  select(-read_count, -gene_id_ms, -gene_len_ms)%>%
  left_join(

readRDS("k19mf/ds/vm-00-sample_info.rds")%>%
  select(sample_id, mouse_num, tumor_type)
)%>%
  relocate(gene_name_ms, gene_name_hu, sample_id, rpkm, mouse_num, tumor_type)




wb <- createWorkbook()

# # Add two sheets to the workbook
addWorksheet(wb, "rpkms")
addWorksheet(wb, "stats")

# Write data frames to sheets
writeData(wb, "rpkms", rpkms)
writeData(wb, "stats", stats)

# Save workbook to an Excel file
saveWorkbook(wb, "k19mf/ds/immune2-out.xlsx", overwrite = TRUE)