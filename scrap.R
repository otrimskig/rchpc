




































plot_list_v1[[1]]


plot_list_v2<-list()



#create grob for thing###############################################
# for (plot_split_index in 1:length(split_dfs)){
#   
#   spec_df1 <- as_tibble(split_dfs[[plot_split_index]]) %>%
#     setNames(gsub("^d\\$", "", names(as_tibble(split_dfs[[plot_split_index]]))))
#   
#   spec_df_filtered <- spec_df1 %>%
#     group_by(across(all_of(split_curves_by))) %>%
#     filter(n() >= 2) %>%
#     ungroup()
#   
#   results <- pairwise_survdiff(
#     reformulate(split_curves_by, response = "Surv(time = aod, event = event)"), 
#     data = spec_df_filtered, 
#     p.adjust.method = "none"
#   )
#   
#   pvals_table<-results[["p.value"]]
#   
#   if(dim(pvals_table)[1]==0){
#     
#     print("no comparisons in table.")
#     
#   }else{

for (plot_split_index in 1:length(split_dfs)) {
  
  spec_df1 <- as_tibble(split_dfs[[plot_split_index]]) %>%
    setNames(gsub("^d\\$", "", names(as_tibble(split_dfs[[plot_split_index]]))))
  
  spec_df_filtered <- spec_df1 %>%
    group_by(across(all_of(split_curves_by))) %>%
    filter(n() >= 2, any(event == 1)) %>%  # Ensure at least one event in each group
    ungroup()
  
  
  
  if (n_distinct(spec_df_filtered[[split_curves_by]]) >= 2) {
    
    results <- pairwise_survdiff(
      reformulate(split_curves_by, response = "Surv(time = aod, event = event)"), 
      data = spec_df_filtered, 
      p.adjust.method = "none"
    )
    
    pvals_table <- results[["p.value"]]
    
    
    if (dim(pvals_table)[1] > 0) {
      
      comps<-results[["p.value"]]%>%
        as_tibble()%>%
        mutate(group_a = attributes(results[["p.value"]])[["dimnames"]][[1]])%>%
        relocate(group_a)%>%
        arrange(group_a)%>%
        pivot_longer(cols=2:last_col(), names_to = "group_b", values_to = "p_value")%>%
        
        mutate(padj_method=results[["p.adjust.method"]])%>%
        
        mutate(grouping_category=results[["grouping"]])%>%
        
        as_tibble()%>%
        filter(!is.na(p_value))%>%
        filter(p_value<=.05)
      
      
      comp_plot<-comps %>%
        filter(p_value<=.05)%>%
        mutate(group_a = factor(group_a, levels = unique(group_a)),
               group_b = factor(group_b, levels = unique(group_b)))%>%
        mutate(row_id = row_number())%>%
        
        
        # Create the comp "plot"
        ggplot(.) +
        geom_tile(aes(x = 7.5, y = 7.5), fill = "#ccffff", alpha=0, width = 15, height = 15)+  # Empty tile
        geom_tile(aes(x = .75, y = row_id, fill = group_a), width = 1.5, height = .85, color="#5c5c5c")+
        geom_tile(aes(x = 2.35, y = row_id, fill = group_b), width = 1.5, height = .85, color="#5c5c5c")+
        
        geom_text(aes(x = 3.3, y = row_id, label = format(p_value, scientific = TRUE, digits = 3), #set rounding
                      
                      fontface = ifelse(p_value <= .05, "bold", "plain"),
                      
        ),
        size=2,
        hjust = 0) +
        
        scale_fill_manual(values = col_map[[split_curves_by]]) +
        theme_minimal() +
        
        theme(
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none",
          plot.margin = margin(0, 0, 0, 0))+
        coord_fixed(ratio = 1)
      
      
    }
  }
  
}

}else{
  print("nrow of df2 is 0.")
}

}else{
  print("subset is na")
}























































plot_list_v2[[plot_split_index]]<-plot_list_v1[[plot_split_index]] +
  annotation_custom(
    grob = ggplotGrob(comp_plot),  # The grob of the second plot
    xmin = -10, xmax = 40,  # Set the x-range for placing the plot inside
    ymin = 0, ymax = .5   # Set the y-range for placing the plot inside
  )


}
}
# 
# if(length(plot_list_v2)>=0){}

plot_list_v3<-list()
#add in meta src data below plots.####################################
for(i in 1:length(plot_list_v2)){
  
  metadata_text <- paste0("src: ", 
                          rstudioapi::getSourceEditorContext()$path %>%
                            sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.) %>%
                            sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.), 
                          "\n", "time: ", 
                          lubridate::round_date(Sys.time(), "second"))
  
  text_grob <- textGrob(metadata_text, x=.05, just="left", gp = gpar(fontsize = 7, col = "gray30"))
  
  plot_list_v3[[i]] <- grid.arrange(
    plot_list_v2[[i]], 
    text_grob, 
    ncol = 1, 
    heights = c(3, 0.3),  # Maintain spacing
    layout_matrix = rbind(c(1), c(2))  # Keeps layout structure stable
  )
  
}


#grid.draw(plot_list_v3[[1]])

##save all plots to pdfs############################################
for(i in 1:length(plot_list_v3)){
  
  ggsave(paste0(plot_output_loc,
                plot_output_base,
                split_plots_by, "-", split_curves_by, "-", plot_subset_clean, "-plot_", subset, "-", i,".pdf"),
         
         title=paste0(gsub("\n", "][", plot_list_v1[[1]]$labels$title), "][",
                      
                      "src: ",
                      rstudioapi::getSourceEditorContext()$path%>%
                        sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                        sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                      
                      " at ", 
                      
                      lubridate::round_date(Sys.time(), "second")
         ),
         
         plot=plot_list_v3[[i]],
         limitsize = FALSE,
         
         height=5,
         width=10,
         scale = 1.2,
         dpi=600,
  )
}

}else{
  print("no data in table. going to next")
}

}else{
  print("na subset.")
}

}

