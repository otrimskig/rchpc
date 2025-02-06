results <- df1 %>%
  group_split(!!sym(split_by))%>%
  
  setNames(sort(unique(df1[[split_by]])))%>%
  
  lapply(function(sub_df, grade) {
    result <- pairwise_survdiff(Surv(time = aod, event = event) ~ patho_cat_name, 
                                data = sub_df, 
                                p.adjust.method = "none")
    
    # Add a new element at the same level as p.value
    result$grouping <- split_by  # This will add 'group' at the same level as other components like p.value
    return(result)
  } )


comps<-list()

for(i in 1:length(names(results))){
  
  category<-paste0(names(results)[i]) 
  
  assign(category, results[[names(results)[i]]])
  
  
  pvals_table<-get(sym(category))[["p.value"]]
  
  if(dim(pvals_table)[1]==0){
    
    print("no comparisons in table.")
    
  }else{
    
    
    #square version with duplicated info (my personal preference)
    comps[[category]]<-get(sym(category))[["p.value"]]%>%
      as_tibble()%>%
      mutate(group_a = attributes(get(sym(category))[["p.value"]])[["dimnames"]][[1]])%>%
      relocate(group_a)%>%
      arrange(group_a)%>%
      pivot_longer(cols=2:last_col(), names_to = "group_b", values_to = "p_value")%>%
      
      mutate(padj_method=get(sym(category))[["p.adjust.method"]])%>%
      
      mutate(grouping_category=get(sym(category))[["grouping"]])%>%
      mutate(grouping_value=sub("cat_","",sym(category)))%>%
      as_tibble()
    
    
  }
  
}

