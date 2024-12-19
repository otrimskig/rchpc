library(tidyverse)
library(ggplot2)
library(ggrepel)
library(edgeR)

r_folder_name<-"acral_paired"
vm02<-readRDS("acral_paired/ds/xv02-filtered_rpkms.rds")

# sample_info<-readRDS("ds/v07-per_sample_info.rds")


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
saveRDS(p, paste0(r_folder_name, "/ds/xmds_dge.rds"))



#create plots folder. 
dir.create(paste0(r_folder_name, "/plots"), showWarnings = FALSE)


pdf(paste0(r_folder_name, "/plots/","xMDS_plot1.pdf"), width = 7, height = 7)  # Set dimensions as needed
plotMDS.DGEList(dge, main = "MDS Plot")
dev.off()



