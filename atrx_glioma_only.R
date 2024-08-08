#this code allows for easy generation of a diff exp dataset between 2 groups. 
#use the generated dataset elsewhere for visualisation. 
#choose the groups and the category and source to generate.

source("libs.R")
library(tidyverse)
library(edgeR)
library(dtplyr)

sample_info<-readRDS("ds/v10-per_sample_updated.rds")

#filter here##


sample_info_fil<-sample_info%>%
  filter(patho_cat=="3")%>%
  filter(resultant_geno=="nf1 KO; pten KO; ink KO; atrx KO"|resultant_geno=="nf1 KO; pten KO; ink KO; atrx wt")



read_counts<-readRDS("ds/vm-02-filtered_rpkms.rds")

all_data<-read_counts%>%
  semi_join(sample_info_fil, by="sample_id")%>%left_join(sample_info_fil)%>%ungroup()%>%arrange(sample_id)

#set category.
#subset.

category<-"resultant_geno"


comp_el<-all_data%>%ungroup()%>%select(!!sym(category))%>%
  count(!!sym(category))%>%select(1)%>%
  pull()


# Assuming comp_el is a factor
if ("NED" %in% levels(comp_el)) {
  # Move "NED" to the first position
  comp_el <- factor(comp_el, levels = c("NED", levels(comp_el)[levels(comp_el) != "NED"]))
}

if ("No evidence of disease" %in% levels(comp_el)) {
  # Move "No evidence of disease" to the first position
  comp_el <- factor(comp_el, levels = c("No evidence of disease", levels(comp_el)[levels(comp_el) != "No evidence of disease"]))
}

if ("nf1 wt; pten wt; ink wt; atrx wt" %in% levels(comp_el)) {
  # Move "nf1 wt; pten wt; ink wt; atrx wt" to the first position
  comp_el <- factor(comp_el, levels = c("nf1 wt; pten wt; ink wt; atrx wt", levels(comp_el)[levels(comp_el) != "nf1 wt; pten wt; ink wt; atrx wt"]))
}



comb_mat<-combn(comp_el, 2)



########################################################################

for(c in 1:ncol(comb_mat)){  

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
        
        
        #save output file.
        saveRDS(output2, paste0("m-dexps/", "dexp-", fs::path_sanitize(paste0("glioma-atrx", category, "-", ga,  " v. ", gb)), ".rds"))
        
        
        
      } else {
        print("issue with group b size.")
      }
      
    } else {  
      print("issue with group a size.")
    }  
  })  
  
  
  
}




