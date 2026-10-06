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


exp_loc<-"/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/exp_data/genewiz/21403"

ref <- readDNAStringSet("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/exp_data/genewiz/refs/nf1crisprfasta_new.fasta")


guide <- DNAString("AAGTCAGCACCGAGCACAAC")


bams<-list.files(paste0(exp_loc,"/results/03_calmd"),
           pattern="filtered.bam$",
           full.names = TRUE)

#x<-20
window<-37

target_start <- 1608 - 37 #window
target_end   <- 1608 + 37 #window

target <- GRanges(
  seqnames = names(ref)[1],
  ranges = IRanges(
    start = target_start,
    end = target_end),
  strand = "-"
  
)

# 
# names(ref)
# width(ref)
cat("\n\n\n\n")

for(bam in 1:length(bams)){
  reads0<-bams[bam]


sample_name<-as.character(sub("_raw\\.calmd\\.bam$","",basename(bams[bam])) )
  
cat("working on", sample_name, "\n\n\n\n")


# 
# matches <- vmatchPattern(
#   guide,
#   ref
# )
# 
# matches


reference <- subseq(
  ref[[1]],
  start = target_start,
  end = target_end
)

reference
#
seqnames(target)

# # Search forward strand
 match_fwd <- matchPattern(
   guide,
   reference
 )
#
# # Search reverse-complement strand
guide_rc <- reverseComplement(guide)
#
match_rev <- matchPattern(
   guide_rc,
   reference
 )
#
match_fwd
match_rev





target0<-target
reference0<-reference
target.loc0<-window-2
  
crispr_set <- readsToTarget(reads0, 
                            target = target0, 
                            reference = reference0,
                            target.loc = target.loc0,
                            #allow.nonnarrowing=TRUE,
                            collapse.pairs = TRUE,
                            #orientation = "positive",
                            orientation = "opposite",
                            chimera.to.target = 200)


# 1. Open the PDF graphics device and specify the filename and dimensions
pdf(paste0("nf1crisp/crispr_variants_plot-",sample_name,".pdf"), width = 18, height = 8)

# 2. Call your plotVariants function
plotVariants(crispr_set)
# 3. Close the device to save the file
dev.off()





cat("finished", sample_name, "\n\n\n\n")



}

# length(reference0)
# target.loc0
# length(target0)
# 
# 
# alns_bad <- readGAlignments(reads0)
# 
# df0<-as.data.frame(alns_bad)
# 
# 
# 
# aln_df <- as.data.frame(alns_bad)
# aln_df[aln_df$width > 900 | grepl("S", aln_df$cigar),
#        c("seqnames", "strand", "cigar", "qwidth", "start", "end", "width")]
# sum(aln_df$width > 900 | grepl("S", aln_df$cigar))
# 
# 
# alns_bad
# alns_bad[2185]
# alns_bad[2186]
# alns_bad[2187]
# 
# 
# aln_df%>%
#   arrange(desc(qwidth))%>%
#   head(30)
# 
# 


install.packages("gridGraphics")

library(gridGraphics)





plotAlignments(crispr_set)
plotFreqHeatmap(crispr_set)
a<-plotVariants(crispr_set)
grid.grab()
a


p_grab <- grid.grab()

grid.draw(p_grab)

library(gridExtra)

title <- textGrob(
  "27634_raw",
  gp = gpar(fontsize = 16, fontface = "bold")
)

info <- textGrob(
  "Guide: AAGTCAGCACCGAGCACAAC | Target: 1608 | Reference: nf1crisprfasta",
  gp = gpar(fontsize = 10)
)

final <- arrangeGrob(
  title,
  info,
  p_grab,
  ncol = 1
)

grid.newpage()
grid.draw(final)





final <- arrangeGrob(
  title,
  info,
  p_grab,
  ncol = 1,
  heights = c(
    0.08,
    0.05,
    0.87
  )
)

grid.newpage()
grid.draw(final)
