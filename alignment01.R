source("libs.R")

library(dplyr)
library(readr)
library(ggplot2)
library(Rsamtools)
library(GenomicAlignments)
library(Biostrings)
library(dplyr)
library(readr)
library(Gviz)
library(jpeg)
library(pwalign)

# Reinstall the CRAN jpeg dependency
#install.packages("jpeg", type = "source")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
#BiocManager::install("pwalign", lib="/uufs/chpc.utah.edu/common/home/u1413890/R/manual_libs")


exp_loc<-"/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/exp_data/genewiz/83338"


variants <- read_tsv(
  paste0(exp_loc,"/results/07_read_analysis/27634_raw.read_variants.tsv"),
  show_col_types = FALSE
)

variants

variants%>%
  filter(ref_pos>1500&ref_pos<1750)%>%
  ggplot(aes(x=ref_pos))+
  geom_bar()




bam <- paste0(exp_loc,"/results/03_calmd/27634_raw.calmd.bam")


ref <- readDNAStringSet(paste0(exp_loc,"/nf1crisprfasta_new.fasta"))
ref_seq <- ref[[1]]





param <- ScanBamParam(
  flag = scanBamFlag(
    isUnmappedQuery = FALSE,
    isSecondaryAlignment = FALSE,
    isSupplementaryAlignment = FALSE
  ),
  mapqFilter = 20,
  what = c(
    "qname",
    "flag",
    "rname",
    "pos",
    "mapq",
    "cigar",
    "seq",
    "qual"
  )
)

bam_data <- scanBam(
  bam,
  param = param
)[[1]]


reads <- tibble(
  read_id  = bam_data$qname,
  flag     = as.integer(bam_data$flag),
  chrom    = bam_data$rname,
  start    = as.integer(bam_data$pos),
  mapq     = as.integer(bam_data$mapq),
  cigar    = bam_data$cigar,
  sequence = as.character(bam_data$seq),
  quality  = as.character(bam_data$qual)
) %>%
  mutate(
    strand = ifelse(
      bitwAnd(flag, 16L) != 0L,
      "-",
      "+"
    )
  )




reads_track <- AlignmentsTrack(
  bam,
  isPaired = FALSE,
  name = "Reads",
  type = "pileup"
)

plotTracks(
  list(
    GenomeAxisTrack(),
    reads_track
  ),
  chromosome = unique(as.character(reads$chrom))[1],
  from = 1,
  to = length(ref_seq)
)


plotTracks(
  list(
    GenomeAxisTrack(),
    reads_track
  ),
  chromosome = "nf1crisprfasta",
  from = 1500,
  to = 1750
)


plotTracks(
  list(
    GenomeAxisTrack(),
    reads_track
  ),
  chromosome = "nf1crisprfasta",
  from = 1600,
  to = 1620,
  main = "Amplicon positions 1500–1750"
)


reads2 <- reads %>%
  mutate(
    ref_span = sapply(cigar, function(x) {
      ops <- str_extract_all(x, "[0-9]+[MDN=X]")[[1]]
      sum(as.integer(str_extract(ops, "[0-9]+")))
    }),
    ref_end = start + ref_span - 1
  )


locus_start <- 1583
locus_end   <- 1633

locus_reads <- reads2 %>%
  filter(
    chrom == "nf1crisprfasta"
  )

nrow(locus_reads)

gal <- readGAlignments(bam, use.names = TRUE)

gal
gal_locus <- subsetByOverlaps(
  gal,
  GRanges(
    seqnames = "nf1crisprfasta",
    ranges = IRanges(1583, 1633)
  )
)

length(gal_locus)











parse_alignment <- function(read_id, start, cigar, sequence, ref_seq, strand) {
  
  # Reverse-complement reverse-strand reads
  seq <- DNAString(sequence)
  
  if (strand == "-") {
    seq <- reverseComplement(seq)
  }
  
  seq <- as.character(seq)
  
  # Parse CIGAR
  ops <- regmatches(
    cigar,
    gregexpr("[0-9]+[MIDNSHP=X]", cigar)
  )[[1]]
  
  cigar_df <- tibble(
    n = as.integer(sub("[A-Z=]", "", ops)),
    op = sub("[0-9]+", "", ops)
  )
  
  ref_pos <- start
  query_pos <- 1
  
  result <- list()
  
  for (i in seq_len(nrow(cigar_df))) {
    
    n <- cigar_df$n[i]
    op <- cigar_df$op[i]
    
    if (op %in% c("M", "=", "X")) {
      
      rp <- ref_pos:(ref_pos + n - 1)
      qp <- query_pos:(query_pos + n - 1)
      
      result[[length(result) + 1]] <- tibble(
        read_id = read_id,
        ref_pos = rp,
        read_base = substr(seq, qp, qp)
      )
      
      ref_pos <- ref_pos + n
      query_pos <- query_pos + n
      
    } else if (op == "I") {
      
      # Put insertion immediately before the next reference base
      qp <- query_pos:(query_pos + n - 1)
      
      result[[length(result) + 1]] <- tibble(
        read_id = read_id,
        ref_pos = ref_pos - 0.1,
        read_base = paste0(
          substr(seq, qp, qp),
          collapse = ""
        )
      )
      
      query_pos <- query_pos + n
      
    } else if (op == "D") {
      
      ref_pos <- ref_pos + n
      
    } else if (op == "S") {
      
      query_pos <- query_pos + n
      
    } else if (op == "N") {
      
      ref_pos <- ref_pos + n
    }
  }
  
  bind_rows(result)
}



locus <- GRanges(
  "nf1crisprfasta",
  IRanges(1583, 1633)
)

param <- ScanBamParam(
  which = locus,
  flag = scanBamFlag(
    isUnmappedQuery = FALSE,
    isSecondaryAlignment = FALSE,
    isSupplementaryAlignment = FALSE
  ),
  mapqFilter = 20
)

stacked <- stackStringsFromBam(
  bam,
  param = param,
  use.names = TRUE
)






stacked


length(stacked)
width(stacked)
head(names(stacked))

as.character(stacked[1:10])




ref_window <- subseq(
  ref_seq,
  start = 1583,
  end = 1633
)

ref_window

data.frame(
  position = 1583:1633,
  reference = as.character(strsplit(as.character(ref_window), "")[[1]]),
  af
)


length(stacked)
width(stacked)
head(names(stacked))
as.character(stacked[1:10])

af <- alphabetFrequencyFromBam(
  bam,
  param = param,
  baseOnly = TRUE
)

af[which(1583:1633 == 1608), ]





library(Biostrings)

ref_window <- subseq(
  ref_seq,
  start = 1583,
  end = 1633
)

ref_set <- DNAStringSet(
  rep(
    as.character(ref_window),
    length(stacked)
  )
)

comparison <- pwalign::compareStrings(
  stacked,
  ref_set
)


head(comparison)



comparison_chr <- as.character(comparison)

read_summary <- data.frame(
  read_id = names(stacked),
  comparison = comparison_chr,
  stringsAsFactors = FALSE
)

read_summary$n_mismatch <- lengths(
  gregexpr("\\?", comparison_chr)
)

read_summary$n_deletion <- lengths(
  gregexpr("-", comparison_chr)
)







c(
  no_difference = sum(read_summary$n_difference == 0),
  altered = sum(read_summary$n_difference > 0)
)




read_summary[
  read_summary$n_difference > 0,
  c("read_id", "comparison", "n_mismatch", "n_deletion")
]

























