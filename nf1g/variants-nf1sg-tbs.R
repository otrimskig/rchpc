
# 
mpileup_file <- sam_puts[23]
cut_site <- 79275040  # Define the CRISPR cut site position
window <- 10          # Define the window around the cut site

mpileup_data <- read.table(mpileup_file, header = FALSE, stringsAsFactors = FALSE)
colnames(mpileup_data) <- c("chrom", "pos", "ref", "depth", "read_bases", "base_qual")

# # Filter for region around cut site
# target_region <- mpileup_data %>%
#     filter(pos >= (cut_site - window) & pos <= (cut_site + window))
#   
#   # Initialize counters
#   insertions <- 0
#   deletions <- 0
#   snps <- 0
#   total_reads <- sum(target_region$depth)       # Count of rows (positions) covered in the target region
#   reads_with_mutations <- 0      
#   reads_with_no_mutations <- 0 # Count of reads with any mutation
#   
#   # Loop over rows and check for SNPs and indels
#   for (i in 1:nrow(target_region)) {
#     i<-1
#     read_bases <- target_region$read_bases[i]
#     ref_base <- target_region$ref[i]
#     depth <- target_region$depth[i]    
#     
#     
#     
#     
#     str_extract_all(read_bases, "\\+\\d+")[[1]]
#     
#     
#     
#     # Count insertions and deletions
#     insertion_matches <- str_extract_all(read_bases, "\\+\\d+")[[1]]
#     deletion_matches <- str_extract_all(read_bases, "-\\d+")[[1]]
#     has_insertion <- length(insertion_matches) > 0
#     has_deletion <- length(deletion_matches) > 0
#     
#     # Count SNPs by comparing non-indel characters to the reference base
#     cleaned_bases <- gsub("\\+\\d+|\\-\\d+", "", read_bases)
#     snp_bases <- str_split(cleaned_bases, "")[[1]]
#     snp_count <- sum(snp_bases %in% c("A", "C", "G", "T") & snp_bases != ref_base)
#     has_snp <- snp_count > 0
#     #no_snp<-snp_count=0
#     
#     
#     # Update individual mutation counts
#     insertions <- insertions + length(insertion_matches)
#     deletions <- deletions + length(deletion_matches)
#     snps <- snps + snp_count
#     
#     # Increment reads_with_mutations if any mutation (insertion, deletion, or SNP) is present
#     if (has_insertion || has_deletion || has_snp) {
#       reads_with_mutations <- reads_with_mutations + 1
#     }
#     
#     if (!has_insertion & !has_deletion & !has_snp) {
#       reads_with_no_mutations <- reads_with_no_mutations + 1
#     }
#     
#     
#   }
#   
#   # Calculate ratio of reads with mutations to those without mutations
#   #reads_without_mutations <- total_reads - reads_with_mutations
#   #mutation_ratio <- ifelse(reads_without_mutations > 0, reads_with_mutations / reads_without_mutations, NA)
#   mutation_perc <- reads_with_mutations/(reads_with_mutations+reads_with_no_mutations)*100
#   
#   # Return results as a list
#   list(
#     cut_site = cut_site,
#     insertions = insertions,
#     deletions = deletions,
#     snps = snps,
#     total_reads = total_reads,
#     reads_with_mutations = reads_with_mutations,
#     #reads_without_mutations = reads_without_mutations,
#     reads_with_no_mutations = reads_with_no_mutations,
#     #mutation_ratio = mutation_ratio,
#     mutation_perc=mutation_perc
#   )
# 
# 
# 
# 











