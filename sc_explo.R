source("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/rchpc/libs.R")


library(tidyverse)
library(Seurat)
library(dplyr)
library(Seurat)
library(purrr)
library(DropletUtils)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(scater)
library(Matrix)
library(DoubletFinder)


setwd("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/exp_data/GSE185386/sc/")

sample_files<-list.files(pattern="MBM01")

pat<-"scRNAseq"
doublet_rate=.049

# 
# 
# c1<-read_csv("GSM6022214_MBM01_sc_all_contig_annotations.csv")
# c2<-read_csv("GSM6022214_MBM01_sc_clonotypes.csv")
# c3<-read_csv("GSM6022214_MBM01_sc_consensus_annotations.csv")
# c4<-read_csv("GSM6022214_MBM01_sc_filtered_contig_annotations.csv")
# 



# Read the data
counts0 <- c10<-read_csv("GSM6022251_MBM01_sc_counts.csv")%>%
  janitor::clean_names()

# Convert to matrix format (gene names as rownames)
gene_names <- counts0[[1]]  # First column is gene names
counts <- counts0[,-1]      # Remove gene name column



counts<-counts0%>%
  column_to_rownames("x1")
  

# Convert to matrix
counts_matrix <- as.matrix.data.frame(counts)





### Identify doublets using scrublet
#doublet_rate_tmp <- doublet_rate[doublet_rate$sample == pat, 2]

doublet_rate_tmp<-doublet_rate



Matrix::writeMM(seu_raw@assays$RNA@counts, paste0('data/', pat, '/matrix_', pat, '_raw.mtx'))


system(paste('python3 brain_mets/scrublet/scrublet_code.py', pat, doublet_rate_tmp))

doublets <- read.table(paste0('data/', pat, '/doublets_', pat, '_raw.txt'), header = T)

seu_raw[['predicted_doublets']] <- doublets$predicted_doublets

seu_raw[['doublet_scores']] <- doublets$doublet_scores
system(paste0('rm data/', pat, '/matrix_', pat, '_raw.mtx'))
system(paste0('rm data/', pat, '/doublets_', pat, '_raw.txt'))





seu_raw <- CreateSeuratObject(counts = counts_matrix, 
                             project = pat, 
                              min.cells = 1, 
                              min.features = 1)



seu_raw <- NormalizeData(seu_raw)
seu_raw <- FindVariableFeatures(seu_raw)
seu_raw <- ScaleData(seu_raw)
seu_raw <- RunPCA(seu_raw)
seu_raw <- RunUMAP(seu_raw, dims = 1:40)
seu_raw <- FindNeighbors(seu_raw, dims = 1:40)
seu_raw <- FindClusters(seu_raw, resolution = 0.1)



### Identify doublets using DF
sweep.list <- paramSweep(seu_raw, PCs = 1:40, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK <- bcmvn %>%
  arrange(desc(BCmetric))
pK <- pK[1, 2]
pK <- as.numeric(levels(pK[[1]]))[pK[[1]]]
nExp <- round(doublet_rate_tmp * dim(seu_raw@assays$RNA@counts)[2])
seu_raw <- doubletFinder_v3(seu_raw, PCs = 1:40, pK = pK, nExp = nExp)
seu_raw$doublet <- seu_raw@meta.data[, paste0('DF.classifications_0.25_', pK, '_', nExp)]
seu_raw$DF_score <- seu_raw@meta.data[, paste0('pANN_0.25_', pK, '_', nExp)]

### Filtering steps 
minFeature <- 500
maxFeature <- 10000
minCount <- 1000
maxCount <- 60000
maxMT <- 10

seu <- subset(seu_raw, subset = nFeature_RNA > minFeature & nFeature_RNA < maxFeature &
                nCount_RNA > minCount & nCount_RNA < maxCount & percent.mt < maxMT & 
                doublet == 'Singlet' & predicted_doublets == F)



DimPlot(seu_raw, reduction = "umap", label = TRUE) +
  ggtitle("UMAP of scRNA-seq data")




# umap_coords <- Embeddings(seurat_obj, "umap")
# write_csv(as.data.frame(umap_coords), "umap_coordinates.csv")



