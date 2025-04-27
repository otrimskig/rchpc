source("libs.R")


library(edgeR)
library(dtplyr)

sample_info<-readRDS("nf1g/ds/v10-per_sample_updated.rds")
read_counts<-readRDS("nf1g/ds/vm-02-filtered_rpkms.rds")


gens<-sample_info$resultant_geno%>%unique()



all_data<-read_counts%>%left_join(sample_info)%>%
  ungroup()%>%arrange(sample_id)%>%
  filter(resultant_geno==gens[2]|resultant_geno==gens[1])





subsetf1<-"patho_cat_name"
subsetf2<-"resultant_geno"


f1s<-all_data%>%
  select(sym(!!subsetf1))%>%unique()%>%
  mutate_if(is.factor, ~droplevels(.x))%>%
  pull()
  





for(f1 in 1:length(f1s)){
  
  cat("\n","working on ", cyan(f1s[f1]), "\n")
  
  
  f1_df0<-all_data%>%
    filter(!!sym(subsetf1)==as.character(f1s[f1]))
  
  


f2s<-f1_df0%>%
  select(sym(!!subsetf2))%>%unique()%>%
  mutate_if(is.factor, ~droplevels(.x))%>%
  pull()


if(length(f2s)<2){
  cat(red("no comparisons to make.\n"))
}else{


      comb_mat <- combn(f2s, 2)
      # Proceed with further processing

      for(c in 1:ncol(comb_mat)){
        ga<-comb_mat[1,c]
        gb<-comb_mat[2,c]
        
        
        #sym() removes quotes. 
        ct<-sym(subsetf2)
        
        
        # !! before category allows it to be evaluated as an object for variable (column) selection. 
        comp_info<-f1_df0%>%
          group_by(sample_id)%>%slice(1)%>%ungroup()%>%
          select(sample_id:last_col())%>%
          select(-read_count, -rpkm, -total_reads)%>%
          
          
          filter(!!ct==ga|!!ct==gb)
        
        group_a_ids<-comp_info%>%
          filter(!!ct==ga)%>%pull(mouse_num)
        
        group_a_count<-comp_info%>%
          filter(!!ct==ga)%>%count()%>%pull()
        
        group_b_ids<-comp_info%>%
          filter(!!ct==gb)%>%pull(mouse_num)
        
        group_b_count<-comp_info%>%
          filter(!!ct==gb)%>%
          count()%>%pull()
        
        tryCatch({
          if (group_a_count>=2) {
            
            cat("group a size checked - ", green("moving on.\n"))
            
            if (group_b_count>=2){
              cat("group b size checked - ", green("moving on.\n"))
              
              cat(blue("reformating data...\n"))
              
              
              comp_counts<-f1_df0%>%
                
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
              
              
              cat("saving output for diff exp...")
              
              
              
              f1_sub_name<-as.character(f1s[f1])
              
              output5<-list(meta=list(src=rstudioapi::getSourceEditorContext()$path%>%
                              sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                              sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                            
                            datetime=lubridate::round_date(Sys.time(), "second"),
                
                            subset1_name=f1_sub_name,
                            
                            comparison=list(groupa=ga, groupb=gb),
                            group_a=group_a_ids,
                            group_a_count=group_a_count,
                            group_b=group_b_ids,
                            group_b_count=group_b_count
                            
                            
                            
                            ),
                            data=output2)
              
              
             
              
              #save output file.
              saveRDS(output5, paste0("nf1g/custom_diff/dexps/", "dexp-", f1s[f1], "-", fs::path_sanitize(paste0(subsetf2, "-", ga,  " v. ", gb)), ".rds"))
              
              cat(bold(green("saved.\n")))
              
            } else {
              cat(red("issue with group b size.\n"))
            }
            
          } else {  
            
            cat(red("issue with group a size.\n"))
            
            
          }  
        })  
        
        
        
      }


}

cat("done.\n")

}
