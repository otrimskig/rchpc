source("libs.R")
library(tidyverse)
library(readr)
library(purrr)
library(DBI)
library(RSQLite)

# read_gmt <- function(file) {
#   lines <- readLines(file)
#   gmt_list <- lapply(lines, function(line) {
#     parts <- strsplit(line, "\t")[[1]]
#     list(term = parts[1], description = parts[2], genes = parts[-(1:2)])
#   })
#   return(gmt_list)
# }
# 
# # Example usage
# gmt_data <- read_gmt("timex/ds/msigdb.v2024.1.Hs.symbols.gmt")
# str(gmt_data)



conn <- dbConnect(RSQLite::SQLite(), "timex/ds/msigdb_v2024.1.Hs.db")
tables<-dbListTables(conn)
db_list<-list()
for (ta in 1:length(tables)){
table_name<-tables[ta]  
  
db_list[[table_name]]<-dbReadTable(conn, table_name)
}


dbDisconnect(conn)

db_list_spec<-db_list

db_list_spec$sqlite_stat1=NULL
db_list_spec$sqlite_stat4=NULL
db_list_spec$species=NULL
db_list_spec$publication_author=NULL
db_list_spec$namespace=NULL
db_list_spec$gene_set_source_member=NULL
db_list_spec$MSigDB=NULL
db_list_spec$gene_set_license=NULL



head(db_list$)%>%view()
