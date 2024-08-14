source("libs.R")
library(tidyverse)

# library(rentrez)






# Get DB (might take few minutes)

library(SRAdb)
sqlfile <- getSRAdbFile(destfile = paste(Sys.Date(), "SRAmetadb.sqlite.gz", sep= '.'))
sra_con<- dbConnect(SQLite(), sqlfile)

# ##  Get exp accession for run accession:
# dbGetQuery(sra_con, "select run_accession, experiment_accession from sra where run_accession = 'SRR1294518'")
#
#
# ## Get more info:
# dbGetQuery(sra_con, "select * from sra where run_accession = 'SRR1294518'")

## Query for more than one run accession id:
dbGetQuery(sra_con, "select run_accession, experiment_accession from sra where run_accession in ('SRR1294518', 'SRR1294519')")


sample_meta<-read_csv("kircher19/by_sample.csv")


sra_ids<-sample_meta%>%
  mutate(sra_id=paste0("'",sra_id, "'"))%>%
  pull(sra_id)%>%
  toString(.)




qu<-dbGetQuery(sra_con, paste0("select run_accession, experiment_accession from sra where run_accession in (",
                           
                           sra_ids,
                           
                           ")"))


saveRDS(qu, "qu.rds")






















# 
# sample_meta<-read_csv("kircher19/by_sample.csv")%>%
#   rename(srx_id=sra_id)
# 
# 
# srx_ids<-sample_meta%>%
#   pull(srx_id)
# 
# 
# 
# 
# 
# 
# 
# 
# print(srx_id)
# 
# 
# # Example SRR ID
# srx_id <- srx_ids[1]
# 
# # Search for the SRX ID in the SRA database
# search_results <- entrez_search(db = "sra", term = srx_id)
# 
# # Fetch the summary to extract the SRR IDs
# summary <- entrez_summary(db = "sra", id = search_results$ids)
# 
# # Extract the SRR IDs from the summary
# srr_ids <- unlist(lapply(summary$Runs, function(run) run$Run@accession))
# 
# # Print the SRR IDs
# print(srr_ids)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Example SRX ID
# srx_id <- "SRX5041829"  # Replace this with your SRX ID
# 
# # Find related SRR IDs using entrez_link
# links <- entrez_link(dbfrom = "sra", id = srx_id, db = "sra")
# 
# # Extract the SRR IDs
# srr_ids <- links$links$sra_sra
# 
# # If there are multiple SRR IDs, fetch their summaries
# if (!is.null(srr_ids) && length(srr_ids) > 0) {
#   srr_summary <- entrez_summary(db = "sra", id = srr_ids)
#   
#   # Extract the Run accessions (SRR IDs)
#   srr_ids_extracted <- sapply(srr_summary, function(x) x$accession)
#   
#   # Print the SRR IDs
#   print(srr_ids_extracted)
# } else {
#   print("No SRR IDs found.")
# }
# 












