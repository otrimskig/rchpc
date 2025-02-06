
source("libs.R")

library(tidyverse)
library(ggpubr)
library(gridExtra)
library(ggkm)
library(scales)

library(survival)
library(survminer)

#reading in current dataset. 
coh1<-readRDS("nf1g/surv/cohorts-2025-01-07.rds")%>%
  mutate(aod = as.numeric(aod))

#final form for survival analysis. 
df0<-coh1%>%
  
  filter(is.na(exclude_from_hist_reason))%>%
  
  
  #keep only necessary vars
  select(exclude, strain_injection, resultant_geno, 
         
          hist,
         
         
         aod, exp_end_date)%>%
  
  mutate(event = if_else(is.na(exp_end_date),1,0))%>%
  select(-exp_end_date)%>%
  
  #making a convenient variable to easily compare nf1 vs. non-nf1s. 
  mutate(nf1=substr(resultant_geno, 1, 6))%>%
  mutate(geno_bg=substr(resultant_geno, 8, nchar(resultant_geno)))%>%
  
  #most conservative exclusion criteria
  filter(is.na(exclude)|exclude>3)%>%
  
  mutate(aod=if_else(event==0, 150,aod))%>%
  
  mutate(hist_cat=substr(hist, 1,1))%>%
  mutate(hist_grade=substr(hist, 3,4))%>%
  mutate(hist_grade=gsub("\\.", "", hist_grade))%>%
  mutate(hist_grade=gsub("15", "1.5", hist_grade))%>%
  mutate(hist_grade=gsub("^d$", "ned", hist_grade))%>%
  
  
  
  mutate(hist_catgrade=paste0(hist_cat, ".", hist_grade))%>%
  mutate(hist_catgrade=if_else(hist_catgrade=="n.ned", "ned", hist_catgrade))%>%
  
  mutate(hist_cat=if_else(hist_cat=="n", "ned", hist_cat))%>%
  mutate(hist_grade=if_else(hist_grade=="N", "no grade assigned", hist_grade))



nf1_renaming_dataset <- readxl::read_excel("nf1g/surv/nf1-renaming dataset.xlsx", 
                                   sheet = "patho_cat")




df1<-df0%>%
  left_join(nf1_renaming_dataset)






aspectratio<-.6


#make empty list
plot_list <- list() 
master_plot_list <- list() 
#overall survival: all cohorts. 




#Plot 01###################################################
#1st element, main plot with all curves. 

split_by<-"hist_grade"

results <- df1 %>%
  group_split(!!sym(split_by)) %>%
  
  setNames(sort(unique(df1[[split_by]])))%>%
  
  lapply(function(sub_df, grade) {
    result <- pairwise_survdiff(Surv(time = aod, event = event) ~ patho_cat_name, 
                                data = sub_df, 
                                p.adjust.method = "none")
    
    # Add a new element at the same level as p.value
    result$grouping <- split_by  # This will add 'group' at the same level as other components like p.value
    return(result)
  } )



results<-results[4]




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




# master_plot_list[[1]]<-ggplot(df1)+
#   geom_km(aes(time = aod, 
#               color=patho_cat_name,
#               status = event),
#           
#           
          
ggplot(df1)+
            geom_km(aes(time = aod, 
                        color=patho_cat_name,
                        status = event),
          
          linewidth=2)+
  
  facet_wrap(vars(hist_grade))+
  
  
  xlim(0,150)+
  
  theme_classic()+
  #theme(legend.position = c(0.2, 0.2))
  
  labs(title="Overall Survival by tumor type, \n separated by Histology grade",
       x = "Days Post Injection",
       y = "% Survival",
       color="Tumor type"
       
  )+
  theme(plot.title = element_text(hjust = 0.5),
        aspect.ratio=aspectratio)









