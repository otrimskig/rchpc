


library(grid)
library(gridExtra)

plotAlignments(crispr_set)

test <- arrangeGrob(
  p,
  ncol = 1
)

grid.newpage()
grid.draw(test)

library("CrispRVariants")
data("gol_clutch1")



p <- plotVariants(gol, plotAlignments.args = list(top.n = 3),
                  plotFreqHeatmap.args = list(top.n = 3),
                  left.plot.margin = ggplot2::unit(c(0.1,0,5,0.2), "lines"))





library("CrispRVariants")
library("rtracklayer")
library("GenomicFeatures")
library("readxl")

# Load the guide location
gd_fname <- system.file(package="CrispRVariants", "extdata/bed/guide.bed")
gd <- rtracklayer::import(gd_fname)
gdl <- resize(gd, width(gd) + 10, fix = "center")

# The saved reference sequence corresponds to the guide 
# plus 5 bases on either side, i.e. gdl
ref_fname <- system.file(package="CrispRVariants", 
                         "extdata/ptena_GRCHz10_ref.rda")
load(ref_fname)

# Load the metadata table, which gives the sample names
md_fname <- system.file(package="CrispRVariants",
                        "extdata/metadata/metadata.xls")
md <- readxl::read_xls(md_fname, 1)

# Get the list of bam files
bam_dir <- system.file(package="CrispRVariants", "extdata/bam")
bam_fnames <- file.path(bam_dir, md$bamfile)

# Check that all files were found
all(file.exists(bam_fnames))








crispr_set <- readsToTarget(bam_fnames, target = gdl, reference = reference,
                            names = md$Short.name, target.loc = 22,
                            verbose = FALSE)










txdb_fname <- system.file("extdata/GRCz10_81_ptena_txdb.sqlite", 
                          package="CrispRVariants")
txdb <- AnnotationDbi::loadDb(txdb_fname)
plotAlignments(crispr_set, txdb = txdb,create.plot = FALSE)

plotFreqHeatmap(crispr_set)

names(plot_data)

do.call(plotAlignments, plot_data)


print(plot_data)
