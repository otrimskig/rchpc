###note: do not do this in parallel. Do in a sequential loop,
###using as many cores as possible (16 for interactive R session)




###be sure to set the it for paired vs single-end reads. 

library(tidyverse)


#set exp path
exp_path<-"../exp_data/kircher19"






# add featurecounts folder if it doesn't exist.
fc_path <- paste0(exp_path, "/featurecounts")

if (!dir.exists(fc_path)) {
  dir.create(fc_path)
  cat("output directory created: ", fc_path)
} else {
  cat("output directory exists: ", fc_path)
}





#gets list of all bam files.
bams<-list.files(paste0(exp_path, "/bams"), pattern = ".bam$", full.names = TRUE)


for(b in 1:length(bams)){
  
  library(tidyverse)
  library(Rsubread)
  
print(paste0(Sys.time(), ": doing fc for ", bams[b]))
  
  
  
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
                       
                      
                       
                       
                      ########
                      #set paired vs. single-end reads
                      
                      # parameters specific to paired end reads
                       isPairedEnd = FALSE,
                       countReadPairs = TRUE,
                       requireBothEndsMapped = FALSE,
                       checkFragLength = FALSE,
                       minFragLength = 50,
                       maxFragLength = 600,
                       countChimericFragments = TRUE,
                       autosort = TRUE,
                       
                       # number of CPU threads
                       nthreads = 16,
                       
                       # read group
                       byReadGroup = FALSE,
                       
                       # report assignment result for each read
                       reportReads = NULL,
                       reportReadsPath = NULL,
                       
                       # miscellaneous
                       # maxMOp = 10,
                       tmpDir = ".",
                       verbose = FALSE)
  
  
  print(paste0(Sys.time(), ": writing output for ", bams[b]))
  

  
#save counts alignments to txt file.
  
  write.table(
    x=data.frame(reads$annotation[,c("GeneID","Length")],
                 reads$counts,
                 stringsAsFactors=FALSE),
    file=paste0(paste0(exp_path, "/featurecounts/"), sub(".bam$", ".FeatureCounts.txt", basename(bams[b]))),
    quote=FALSE,
    sep="\t",
    row.names=FALSE)
  
#save stats for alignment performance. 
  write.table(
    x=data.frame(reads$stat,
                 stringsAsFactors=FALSE),
    file=paste0(paste0(exp_path, "/featurecounts/"), sub(".bam$", ".FeatureCounts_stats.txt", basename(bams[b]))),
    quote=FALSE,
    sep="\t",
    row.names=FALSE)
  
  
}






