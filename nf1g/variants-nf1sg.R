# Load required libraries
library(tidyverse)
library(dplyr)
library(stringr)

library(dplyr)
library(stringr)

# Function to parse mpileup and count per-read mutations within a specified region
parse_mpileup_for_read_mutations <- function(mpileup_file, cut_site, window = 5) {
  # Read in mpileup file
  mpileup_data <- read.table(mpileup_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(mpileup_data) <- c("chrom", "pos", "ref", "depth", "read_bases", "base_qual")
  
  # Filter for region around cut site
  target_region <- mpileup_data %>%
    filter(pos >= (cut_site - window) & pos <= (cut_site + window))
  
  # Initialize counters
  reads_with_insertions <- 0
  reads_with_deletions <- 0
  reads_with_snps <- 0
  reads_without_mutations <- 0
  
  # Loop over each position in the target region
  for (i in 1:nrow(target_region)) {
    read_bases <- target_region$read_bases[i]
    ref_base <- target_region$ref[i]
    
    # Split read bases into individual read symbols for mutation checking
    individual_bases <- strsplit(read_bases, "")[[1]]
    
    # Initialize flag for each read to track if a mutation was found
    mutation_found <- FALSE
    
    j <- 1
    while (j <= length(individual_bases)) {
      base <- individual_bases[j]
      
      # Check if current base represents an insertion (e.g., "+<number><bases>")
      if (base == "+") {
        reads_with_insertions <- reads_with_insertions + 1
        mutation_found <- TRUE
        # Skip over insertion bases
        j <- j + 2 + as.integer(individual_bases[j + 1])
      }
      # Check if current base represents a deletion (e.g., "-<number><bases>")
      else if (base == "-") {
        reads_with_deletions <- reads_with_deletions + 1
        mutation_found <- TRUE
        # Skip over deletion bases
        j <- j + 2 + as.integer(individual_bases[j + 1])
      }
      # Check if current base is an SNP
      else if (base %in% c("A", "C", "G", "T") && base != ref_base) {
        reads_with_snps <- reads_with_snps + 1
        mutation_found <- TRUE
      }
      
      j <- j + 1
    }
    
    # If no mutation was found in this read, count it as a non-mutated read
    if (!mutation_found) {
      reads_without_mutations <- reads_without_mutations + 1
    }
  }
  
  # Return results as a list
  list(
    cut_site = cut_site,
    reads_with_insertions = reads_with_insertions,
    reads_with_deletions = reads_with_deletions,
    reads_with_snps = reads_with_snps,
    reads_without_mutations = reads_without_mutations,
    percent_mutated = sum(reads_with_deletions,reads_with_snps,reads_with_insertions)/sum(reads_with_deletions,reads_with_snps,reads_with_insertions, reads_without_mutations)*100
  )
}










# Define the file paths and cut site information
sam_puts <- list.files("../exp_data/23908R_merged/spec_counts", pattern = "*.sam_pileup.txt$", full.names = TRUE) %>%
  tibble("path" = .) %>%
  filter(!grepl("*.23908X15.*", path)) %>%
  filter(!grepl("*.23908X07.*", path)) %>%
  pull()














read_variant_counts <- list()

# Process each mpileup file
for (sf in seq_along(sam_puts)) {
  mpileup_file <- sam_puts[sf]
  
  # Extract the sample name from the filename
  sample_name <- basename(mpileup_file) %>%
    sub("_.*", "", .)
  
  # Define the CRISPR cut site position and window
  cut_site <- 79275040
  window <- 10
  
  # Get variant counts per read
  variant_counts <- parse_mpileup_for_read_mutations(mpileup_file, cut_site, window)
  
  # Add sample name and window to the result
  variant_counts$sample <- sample_name
  variant_counts$window_bp <- window*2
  read_variant_counts[[sf]] <- variant_counts
}

# Convert the list to a tibble
read_variant_counts_tibble <- bind_rows(read_variant_counts) %>%
  select(sample, cut_site, window_bp, reads_with_insertions, reads_with_deletions, reads_with_snps, reads_without_mutations)

# Display the combined tibble
print(read_variant_counts_tibble)


b<-read_variant_counts_tibble%>%
  mutate(total_reads=reads_with_insertions + reads_with_deletions + reads_with_snps + reads_without_mutations)%>%
  mutate(reads_with_mutations = total_reads-reads_without_mutations)%>%
  relocate(reads_with_mutations, .before = "reads_without_mutations")%>%
  mutate(mutation_perc=reads_with_mutations/total_reads*100)





saveRDS(b, "nf1g/ds/mutations_at_nf1_guide_site.rds")
