library(foreach)

if (!exists("n.cores")) {
  
  "initilizing cores..."
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  doParallel::registerDoParallel(cl = my.cluster)
  
  #check if it is registered (optional)
  foreach::getDoParRegistered()
  
  "parallel cores initialized."
  
}


bams<-list.files("../exp_data/23908R_merged/", pattern = ".bam$", full.names = TRUE)

foreach(b=1:length(bams)) %dopar% {
  
  library(tidyverse)
  library(Rsubread)
  
  
  
  
  
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
                       # maxMOp = 10,
                       tmpDir = ".",
                       verbose = FALSE)
  
  
#save counts alignments to txt file.
  
  write.table(
    x=data.frame(reads$annotation[,c("GeneID","Length")],
                 reads$counts,
                 stringsAsFactors=FALSE),
    file=paste0("../exp_data/23908R_merged/featurecounts/", sub(".bam$", ".FeatureCounts.txt", basename(bams[b]))),
    quote=FALSE,
    sep="\t",
    row.names=FALSE)
  
#save stats for alignment performance. 
  write.table(
    x=data.frame(reads$stat,
                 stringsAsFactors=FALSE),
    file=paste0("../exp_data/23908R_merged/featurecounts/", sub(".bam$", ".FeatureCounts_stats.txt", basename(bams[b]))),
    quote=FALSE,
    sep="\t",
    row.names=FALSE)
  
  
}






