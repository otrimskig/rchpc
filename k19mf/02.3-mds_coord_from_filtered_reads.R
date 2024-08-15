source("libs.R")

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(edgeR)


#set exp directory, relative to rchpc proj.
exp_dir<-"k19mf/"


vm02<-readRDS(paste0(exp_dir, "ds/vm-02-filtered_rpkms.rds"))

mat<-vm02%>%
  select(gene_id_ms, sample_id, rpkm)%>%
  pivot_wider(names_from = sample_id, values_from = rpkm)%>%
  column_to_rownames("gene_id_ms")%>%
  as.matrix.data.frame()

lib_size<-vm02%>%
  select(sample_id, read_count)%>%
  ungroup()%>%
  group_by(sample_id)%>%
  summarise(lib_size=sum(read_count))%>%
  ungroup()%>%
  arrange(sample_id)%>%
  pull(lib_size)

dge<-DGEList(counts = mat, 
        lib.size = lib_size)



p<-plotMDS.DGEList(dge)

plotMDS.DGEList(dge)


pdf(paste0(exp_dir, "plots/mds_plot-edger.pdf"), width = 12, height = 8)

# Create the plot
plotMDS.DGEList(dge)

# Close the PDF device
dev.off()



saveRDS(p, paste0(exp_dir, "ds/mds_dge.rds"))



