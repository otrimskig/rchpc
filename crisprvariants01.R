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
library(CrispRVariants)
# Reinstall the CRAN jpeg dependency
#install.packages("jpeg", type = "source")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
#BiocManager::install("CrispRVariants", lib="/uufs/chpc.utah.edu/common/home/u1413890/R/manual_libs")


exp_loc<-"/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/exp_data/genewiz/83338"

ref <- readDNAStringSet(
  paste0(exp_loc, "/nf1crisprfasta_new.fasta")
)


guide <- DNAString("AAGTCAGCACCGAGCACAAC")


bams<-list.files("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/exp_data/genewiz/83338/results/03_calmd",
           pattern=".bam$",
           full.names = TRUE)




# 
# names(ref)
# width(ref)

for(bam in 1:length(bams)){
  reads0<-bams[bam]

  
sample_name<-as.character(sub("_raw\\.calmd\\.bam$","",basename(bams[bam])) )
  
  


# 
# matches <- vmatchPattern(
#   guide,
#   ref
# )
# 
# matches
target_start <- 1608 - 10
target_end   <- 1608 + 15

target <- GRanges(
  seqnames = names(ref)[1],
  ranges = IRanges(
    start = target_start,
    end = target_end
  )
)

reference <- subseq(
  ref[[1]],
  start = target_start,
  end = target_end
)

# reference
# 
# seqnames(target)

# # Search forward strand
# match_fwd <- matchPattern(
#   guide,
#   reference
# )
# 
# # Search reverse-complement strand
# guide_rc <- reverseComplement(guide)
# 
# match_rev <- matchPattern(
#   guide_rc,
#   reference
# )
# 
# match_fwd
# match_rev





target0<-target
reference0<-reference
target.loc0<-10
  
crispr_set <- readsToTarget(reads0, 
                            target = target0, 
                            reference = reference0,
                            target.loc = target.loc0,
                            collapse.pairs = TRUE,
                            chimera.to.target = 200)







# 1. Open the PDF graphics device and specify the filename and dimensions
pdf(paste0("nf1crisp/crispr_variants_plot-",sample_name,".pdf"), width = 12, height = 8)

# 2. Call your plotVariants function
plotVariants(crispr_set)
# 3. Close the device to save the file
dev.off()

}

