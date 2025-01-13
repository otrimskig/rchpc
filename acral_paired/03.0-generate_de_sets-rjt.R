##########################
source("libs.R")
library(tidyverse)
library(edgeR)
library(dtplyr)

sample_info<-readRDS("acral_paired/ds/v00-sample_info.rds")
read_counts<-readRDS("acral_paired/ds/v02-filtered_rpkms.rds")

all_data<-read_counts%>%left_join(sample_info)%>%ungroup()%>%arrange(sample_id)%>%
  filter(mouse_num!="38536")

#add things to file name. Be sure to add "_" or "-" pre/post as needed.
manual_file_prefix<-NULL
manual_file_suffix<-"_RJT_exclude38536"



#view all available categories.
colnames(all_data)


#use to determine category for selection of comparison.
categories<-c("sample_type")

r_folder_name<-"acral_paired"


#generate de sets from all categories specified above.
tryCatch({
for(i in 1:length(categories)){  
# i<-1

category<-categories[i]


comp_el<-all_data%>%ungroup()%>%select(!!sym(category))%>%
  count(!!sym(category))%>%select(1)%>%
  pull()

comb_mat<-combn(comp_el, 2)



########################################################################




for (c in 1:ncol(comb_mat)){  

# c<-1 

# sample_info<-readRDS("ds/v10-per_sample_updated.rds")
# read_counts<-readRDS("ds/vm-02-filtered_rpkms.rds")
#   
# all_data<-read_counts%>%left_join(sample_info)%>%ungroup()%>%arrange(sample_id)
  
  
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

  
  
  dir.create(paste0(r_folder_name, "/dexps"), showWarnings = F)
  
  #save output file.
  saveRDS(output2, paste0(r_folder_name, "/dexps/", "dexp-", fs::path_sanitize(paste0(manual_file_prefix, category, "-", ga,  " v. ", gb, manual_file_suffix)), ".rds"))
  
  
  
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






