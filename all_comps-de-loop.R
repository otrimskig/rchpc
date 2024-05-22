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


#this code allows for easy generation of a diff exp dataset between 2 groups. 
#use the generated dataset elsewhere for visualisation. 
#choose the groups and the category and source to generate.

library(tidyverse)
library(edgeR)
library(foreach)

#read in all data. tidy dataframe, with all sample info
#plus every raw count for each sample and gene.
#remove 15 and 21 for low quality. 
all_data<-readRDS("ds/v06-all_counts_plus_info.rds")%>%
  arrange(sample_id)


##########################################################################
#view all available categories.
colnames(all_data)


#use to determine category for selection of comparison.
category<-"resultant_geno"

comp_el<-all_data%>%select(!!sym(category))%>%
  count(!!sym(category))%>%select(1)%>%
  pull()


if ("NED" %in% comp_el) {
  # Move "NED" to the first position
  comp_el <- c("NED", comp_el[comp_el != "NED"])
}

comb_mat<-combn(comp_el, 2)



########################################################################

foreach(c=1:ncol(comb_mat)) %dopar% {  
  
  library(tidyverse)
  library(edgeR)

  
  all_data<-readRDS("ds/v06-all_counts_plus_info.rds")%>%
    arrange(sample_id)
  
  
  ##########################################################################


ga<-comb_mat[1,c]
gb<-comb_mat[2,c]


#sym() removes quotes. 
ct<-sym(category)


# !! before category allows it to be evaluated as an object for variable (column) selection. 
comp_info<-all_data%>%
  group_by(sample_id)%>%slice(1)%>%ungroup()%>%
  select(sample_id:last_col())%>%
  select(-raw_read_count, -rpkm)%>%
  
  
  filter(!!ct==ga|!!ct==gb)



comp_counts<-all_data%>%
  
  #use comp_info to subset all reads by samples included in comparison. 
  semi_join(comp_info, by="sample_id")%>%
  select(gene_name, gene_id, raw_read_count, sample_id)%>%
  
  #pivot df to wider format for input into edgeR. 
  pivot_wider(names_from = sample_id,
              values_from = raw_read_count)%>%
  
  
  #ensure gene_ids are unique before setting rownames.
  group_by(gene_name)%>%slice(1)%>%ungroup()%>%
  column_to_rownames("gene_name")


##################################################

#make data formatted for dge. 
c_counts<-comp_counts%>%select(-gene_id)

c_genes<-comp_counts%>%select(gene_id)

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
  rownames_to_column("gene_name")%>%
  as_tibble()
  


output2<-output1%>%
  left_join(all_data%>%select(gene_name, rpkm, sample_id)%>%
              semi_join(comp_info, by="sample_id"))%>%
  pivot_wider(values_from = rpkm, names_from = sample_id, names_prefix = "rpkm_")




#save output file.
saveRDS(output2, paste0("dexps/", "dexp-", category, "-", ga,  "v", gb, ".rds"))


}
