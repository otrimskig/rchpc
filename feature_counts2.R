setwd("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig")


library(tidyverse)
library(Rsubread)







bams<-list.files("23908R/starbam", pattern = ".bam$", full.names = TRUE)

bams[1:2]

for (b in 1:length(bams)){
  
  reads<-featureCounts(bams[b],
                       
                       # annotation
                       annot.inbuilt = "mm39",
                       annot.ext = NULL,
                       isGTFAnnotationFile = FALSE,
                       GTF.featureType = "exon",
                       GTF.attrType = "gene_id",
                       GTF.attrType.extra = NULL,
                       chrAliases = NULL,
                       
                       # level of summarization
                       useMetaFeatures = TRUE,
                       
                       # overlap between reads and features
                       allowMultiOverlap = FALSE,
                       minOverlap = 1,
                       fracOverlap = 0,
                       fracOverlapFeature = 0,
                       largestOverlap = FALSE,
                       nonOverlap = NULL,
                       nonOverlapFeature = NULL,
                       
                       # Read shift, extension and reduction
                       readShiftType = "upstream",
                       readShiftSize = 0,
                       readExtension5 = 0,
                       readExtension3 = 0,
                       read2pos = NULL,
                       
                       # multi-mapping reads
                       countMultiMappingReads = TRUE,
                       
                       # fractional counting
                       fraction = FALSE,
                       
                       # long reads
                       isLongRead = FALSE,
                       
                       # read filtering
                       minMQS = 0,
                       splitOnly = FALSE,
                       nonSplitOnly = FALSE,
                       primaryOnly = FALSE,
                       ignoreDup = FALSE,
                       
                       # strandness
                       strandSpecific = 0,
                       
                       # exon-exon junctions
                       juncCounts = FALSE,
                       genome = NULL,
                       
                       # parameters specific to paired end reads
                       isPairedEnd = TRUE,
                       countReadPairs = TRUE,
                       requireBothEndsMapped = FALSE,
                       checkFragLength = FALSE,
                       minFragLength = 50,
                       maxFragLength = 600,
                       countChimericFragments = TRUE,
                       autosort = TRUE,
                       
                       # number of CPU threads
                       nthreads = 10,
                       
                       # read group
                       byReadGroup = FALSE,
                       
                       # report assignment result for each read
                       reportReads = NULL,
                       reportReadsPath = NULL,
                       
                       # miscellaneous
                       maxMOp = 10,
                       tmpDir = ".",
                       verbose = FALSE)
  
  
  
  
  write.table(
    x=data.frame(reads$annotation[,c("GeneID","Length")],
                 reads$counts,
                 stringsAsFactors=FALSE),
    file=paste0("23908R/featurecounts/", sub(".out.bam$", ".FeatureCounts.txt", basename(bams[b]))),
    quote=FALSE,
    sep="\t",
    row.names=FALSE)
  
  
  
  
  
  
  
}






combined_data <- read.table(counts_txts[1], header = TRUE)



counts_txts<-list.files("23908R/featurecounts", full.names = TRUE)


for (t in 2:length(counts_txts)){
  
  data <- read.table(counts_txts[t], header = TRUE)
  
  combined_data <- full_join(combined_data, data)
  
  
}



write.table(combined_data, file = "23908R/v01-all_counts.txt", sep = "\t", row.names = FALSE)








