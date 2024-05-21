setwd("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/exp_data/")


library(tidyverse)

library(Rsubread)



reads<-featureCounts("23908R1/bams/23908X01_20240412_LH00227_0067_A22KHYJLT3_S18_L003_R.resultsAligned.sortedByCoord.out.bam",

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
              #maxMOp = 10,
              tmpDir = ".",
              verbose = FALSE)


reads$stat


write.table(
  x=data.frame(reads$stat,
               stringsAsFactors=FALSE),
  file=paste0("../exp_data/23908R_merged/featurecounts/", "FeatureCounts_stats.txt"),
  quote=FALSE,
  sep="\t",
  row.names=FALSE)






