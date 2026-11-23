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
library(gridGraphics)
library(gridExtra)

# file locations. 
exp_loc<-"/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/exp_data/genewiz/21403"

sample_info<-read_csv("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/exp_data/genewiz/21403/sample_info.csv")%>%
  mutate(ms_num=as.character(ms_num))

ref <- readDNAStringSet("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/exp_data/genewiz/refs/nf1crisprfasta_new.fasta")

bams<-list.files(paste0(exp_loc,"/results/03_calmd"),
           pattern="filtered.bam$",
           full.names = TRUE)

#========================================
# set inputs. 
#========================================

#x<-20
window<-37

target_start <- 1608 - 37 #window
target_end   <- 1608 + 37 #window

target0 <- GRanges(
  seqnames = names(ref)[1],
  ranges = IRanges(
    start = target_start,
    end = target_end),
  strand = "-"
  
)

target.loc0<-window-2


#for reference only, won't be used to make plot.
guide <- DNAString("AAGTCAGCACCGAGCACAAC")

#=======================================
# make plots.
#=======================================


cat("\n\n\n\n")
for(bam in 1:length(bams)){

reads0<-bams[bam]
#bam<-9
sample_name<-as.character(sub("_raw\\.filtered\\.bam$","",basename(bams[bam])) )
  
cat("working on", sample_name, "\n\n\n\n")

reference0 <- subseq(
  ref[[1]],
  start = target_start,
  end = target_end
)

crispr_set <- readsToTarget(reads0, 
                            target = target0, 
                            reference = reference0,
                            target.loc = target.loc0,
                            #allow.nonnarrowing=TRUE,
                            collapse.pairs = TRUE,
                            #orientation = "positive",
                            orientation = "opposite",
                            chimera.to.target = 200)




plotVariants(crispr_set)
p_grab <- grid.grab()




#==============================================
#==============================================

#grid.draw(p_grab)

sa<-sample_info%>%
  filter(ms_num==sample_name)

ms_num<-sa%>%pull(ms_num)
strain<-sa%>%pull(strain)
injection<-sa%>%pull(inj)




title <- textGrob(
  paste(ms_num, "\n", strain, "\n", injection),
  gp = gpar(fontsize = 16, fontface = "bold")
)

# info <- textGrob(
#   paste(strain, "\n", injection),
#   gp = gpar(fontsize = 10)
# )




final <- arrangeGrob(
  title,
  p_grab,
  ncol = 1,
  heights = c(
    0.1,
    0.9
  )
)





grid.newpage()

# 1. Open the PDF graphics device and specify the filename and dimensions
pdf(paste0("nf1crisp/crispr_variants_plot-",sample_name,".pdf"), width = 18, height = 8)

grid.draw(final)


# 2. Call your plotVariants function
#plotVariants(crispr_set)
# 3. Close the device to save the file
dev.off()



cat("finished", sample_name, "\n\n\n\n")
dev.off()
}





plots<-fs::dir_info("nf1crisp", full.names = T, 
                    pattern=".pdf$",
                    #pattern="^hm.*\\.pdf$", 
                    recurse = F)%>%
  tibble()%>%
  dplyr::filter(modification_time>Sys.time()-lubridate::minutes(5))%>%
  pull(path)

qpdf::pdf_combine(plots, output = "nf1crisp/plots-combined-1.pdf")






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












# reference
# #
# seqnames(target)
# 
# # # Search forward strand
#  match_fwd <- matchPattern(
#    guide,
#    reference
#  )
# #
# # # Search reverse-complement strand
# guide_rc <- reverseComplement(guide)
# #
# match_rev <- matchPattern(
#    guide_rc,
#    reference
#  )
# #
# match_fwd
# match_rev



