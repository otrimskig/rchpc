##########################
source("libs.R")
library(tidyverse)
library(edgeR)
library(dtplyr)
library(ggplot2)

vm0<-readRDS("acral_sub_rppa/ds/rppa_df_list0.rds")
sample_info0<-readRDS("acral_sub_rppa/ds/sample_info0.rds")
sample_interest<-sample_info0%>%
  filter(sample_type!="none")

vm1<-vm0[["L4 (linear)"]][["df_long"]]
vm2<-vm1%>%
  filter(order %in% sample_interest$order)%>%
  mutate(rppa_value=as.numeric(rppa_value))%>%
  
  left_join(sample_info0%>%select(sample_num, order, sample_type))
  
  # filter(order!=455&order!=456)%>%
  # filter(order!=469)%>%
  
 
saveRDS(vm2, "acral_sub_rppa/ds/all_data_long.rds")
  

# library
library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)
  # 
  # left_join(sample_interest%>%select(mouse_num, order))%>%
  # mutate(sample_name_hm=paste0(mouse_num, "-", order))%>%
  # select(-order, -mouse_num)%>%
  # pivot_wider(names_from = antibody_name, values_from=rppa_value)%>%
  # column_to_rownames("sample_name_hm")%>%
  # as.matrix.data.frame()%>%
  # t()
  # 



rep_diff <- vm2 %>%
  group_by(antibody_name, sample_num, sample_type) %>%
  summarise(
    val1 = rppa_value[1],
    val2 = rppa_value[2],
    mean=mean(c(val1, val2)),
    abs_diff = abs(val1 - val2),
    rel_diff = abs(val1 - val2) / mean(c(val1, val2)),
    .groups = "drop"
  )


vm3<-vm2%>%
  left_join(rep_diff)






# make sure sample_num is a factor for nice ordering
vm3$sample_num <- as.factor(vm3$sample_num)

# plot
ggplot(vm3, aes(x = rel_diff, y = sample_num, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis(name = "Relative Diff", option = "C") +
  labs(
    title = "Distribution of Relative Differences by Sample",
    x = "Relative Difference",
    y = "Sample"
  ) +
  theme_ipsum() +
  theme(
    legend.position = "right",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )










ggplot(vm3, aes(x = rel_diff)) +
  geom_density(linewidth=1,alpha=.2) +
  
  theme_minimal()
  # #labs(
  #   title = "Density of Relative Differences by Sample",
  #   x = "Relative Difference",
  #   y = "Density",
  #   color = "Sample Number"
  # )









filtered_diffs<-vm3%>%
  mutate(rppa_filtered=if_else(rel_diff<=.375, rppa_value, NA))%>%
  mutate(mean_filtered=if_else(rel_diff<=.375, mean, NA))





hm_mat<-filtered_diffs%>%
  mutate(sample_num_type=paste0(sample_num, "_", sample_type))%>%
  select(antibody_name, sample_num_type, mean_filtered)%>%
  unique()%>%
  pivot_wider(names_from = sample_num_type, values_from=mean_filtered)%>%
  column_to_rownames("antibody_name")%>%
  as.matrix.data.frame()


n_threshold <- 9  # for example, keep rows with at least 5 non-NA entries

# Apply filtering to the matrix
hm_mat_filtered <- t(scale(t(hm_mat[rowSums(!is.na(hm_mat)) >= n_threshold, ])))



Heatmap(hm_mat_filtered)







fdf2<-filtered_diffs%>%
  select(-rppa_value, -rppa_filtered, -order)%>%
  unique()
  
library(dplyr)
library(DESeq2) # For the exact test








stop()










########################################################################




foreach(c=1:ncol(comb_mat)) %dopar% {  
# for(c in 1:ncol(comb_mat)){  
 
  source("libs.R")
  library(tidyverse)
  library(edgeR)
  library(dtplyr)
  
# c<-1 

sample_info<-readRDS("ds/v10-per_sample_updated.rds")
read_counts<-readRDS("ds/vm-02-filtered_rpkms.rds")
  
all_data<-read_counts%>%left_join(sample_info)%>%ungroup()%>%arrange(sample_id)
  
  
  ##########################################################################
  
  
  ga<-comb_mat[1,c]
  gb<-comb_mat[2,c]
  
# ga
# gb
  
  #sym() removes quotes. 
  ct<-sym(category)
  
  
  # !! before category allows it to be evaluated as an object for variable (column) selection. 
  comp_info<-all_data%>%
    group_by(sample_id)%>%slice(1)%>%ungroup()%>%
    select(sample_id:last_col())%>%
    select(-read_count, -rpkm, -total_reads)%>%
    
    
    filter(!!ct==ga|!!ct==gb)
  

group_a_count<-comp_info%>%
  filter(!!ct==ga)%>%count()%>%pull()

group_b_count<-comp_info%>%
  filter(!!ct==gb)%>%
  count()%>%pull()

tryCatch({
if (group_a_count>=2) {
  
  print("group a size checked - moving on.")
  
  if (group_b_count>=2){
  print("group b size checked - moving on.")




comp_counts<-all_data%>%

    #use sample_ids contained in comp_info to subset all reads to only samples included in comparison. 
    semi_join(comp_info, by="sample_id")%>%

    #select columns needed for comparison.
    select(gene_name_ms, gene_id_ms, read_count, sample_id)%>%
    
    #pivot df to wider format for input into edgeR. 
    pivot_wider(names_from = sample_id,
                values_from = read_count)%>%
    
    filter(!is.na(gene_name_ms))%>%
    #ensure gene_ids are unique before setting rownames.
    group_by(gene_name_ms)%>%slice(1)%>%ungroup()%>%
    group_by(gene_id_ms)%>%slice(1)%>%ungroup()%>%
    column_to_rownames("gene_name_ms")
  
  
  ##################################################
  
  #make data formatted for dge. 
  c_counts<-comp_counts%>%select(-gene_id_ms)
  
  c_genes<-comp_counts%>%select(gene_id_ms)
  
  c_group<-comp_info%>%
    select(all_of(ct))%>%
    mutate(model_bin=if_else(!!ct==ga, 0, 
                             if_else(!!ct==gb, 1, NA_real_)))%>%
    pull(model_bin)
  
  
  #############################################
  #make dge object and run analysis.
  dge<-DGEList(counts=c_counts, 
               genes =c_genes,
               
               group = c_group)
  
  
  #run standard tests and filters.
  keep <- filterByExpr(dge)
  
  #filter out some genes. 
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  
  #calculate normalization for gene counts. 
  dge<-calcNormFactors(dge)
  
  #stats
  dge<-estimateDisp(dge, design = model.matrix(~c_group))
  
  exacDGE<-exactTest(dge, pair = 1:2, dispersion = "auto", rejection.region = "doubletail")
  
  top<-topTags(exacDGE, n=100000)
  
  
  ####################################
  #reformat output of analysis, and rejoin to other info.
  output1<-as.data.frame(top)%>%
    rownames_to_column("gene_name_ms")%>%
    as_tibble()
  
  
  
  output2<-output1%>%
    left_join(all_data%>%select(gene_name_ms, rpkm, sample_id)%>%
                semi_join(comp_info, by="sample_id"))%>%

    pivot_wider(values_from = rpkm, names_from = sample_id, names_prefix = "rpkm_")

  stop()
  #save output file.
  saveRDS(output2, paste0("m-dexps/", "dexp-", fs::path_sanitize(paste0(category, "-", ga,  " v. ", gb)), ".rds"))
  
  
  
  } else {
    print("issue with group b size.")
  }
  
} else {  
  print("issue with group a size.")
}  
})  
  
  

}
}
  
}) 

