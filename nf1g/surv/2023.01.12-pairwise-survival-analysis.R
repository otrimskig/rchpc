####################################
#2023.01.12
#creating function to easily export survival analysis data. 

#2023.03.08 - updated comments to increase usability

library(tidyverse)
library(survminer)
library(survival)
library(googlesheets4)



#############
#set parameters for analysis in the input and output parameters sections. 
#if wanting a different output style, select version by 
#commenting in/out respective output table style in the processing section. 


##########################################################
#input parameters

#import data into R and make it into convenient form for processing.
#:::generic data format:::

#column for group or any other grouping factor
#column for id of individual
#column for age_of_death / time of event
#column for "event" in binary form, 1/0. 


readRDS("ds/cohort3_survival2.rds")->cohort_survival



#gets member count for each group. 
#assigns each group a number for easier processing later. 

#merges that to imported dataset, such that each individual 
#obs also contains information about the group it belongs to. 

cohort_survival%>%
  count(genes_ko)%>%
  mutate(group = 1:6)->genes_ko

cohort_survival%>%
  left_join(genes_ko, by = "genes_ko")->cohort_survival2



#########################################################
#output/analysis parameters

#set results file name

"cohort_surv_analysis" ->name_of_file


#create survival diff object from imported data. 
#####set methods and variables for analysis here. 

pairwise_survdiff(Surv(time = age_death2, 
                  
                  #event information. 
                  event = event)
                 
                  #grouped by
                  ~ genes_ko, 
                  
                  #dataset
                  data = cohort_survival2, 
                  
                  #adjustment method. 
                  #Allowed values include "holm", "hochberg", "hommel", 
                  #"bonferroni", "BH", "BY", "fdr", "none". 
                  p.adjust.method = "none"
                  
                  )->pairwise


############################################################
#############################################################
#processing 

#get names of groups and tests to be added to the results table below. 
attributes(pairwise[["p.value"]])[["dimnames"]][[1]]->group_a



#choose output version:


#3-column version:
# pairwise[["p.value"]]%>%
#   as_tibble()%>%
#   mutate(group_a = group_a)%>%
#   relocate(group_a)%>%
#   arrange(group_a)%>%
#   pivot_longer(cols=2:6, names_to = "group_b", values_to = "p_value")%>%
#   filter(!is.na(p_value))->results



#square comparison version:
# pairwise[["p.value"]]%>%
#   as_tibble()%>%
#   mutate(group = group_a)%>%
#   relocate(group)%>%
#   arrange(group)->results
  
  


#square version with duplicated info (my personal preference)
pairwise[["p.value"]]%>%
  as_tibble()%>%
  mutate(group_a = group_a)%>%
  relocate(group_a)%>%
  arrange(group_a)%>%
  pivot_longer(cols=2:6, names_to = "group_b", values_to = "p_value")%>%
  mutate(group_c = group_a)->table

table%>%
  select(group_a, group_b, p_value)->tableab

table%>%
  select(group_b, group_c, p_value)%>%
  rename(group_a=group_b, group_b=group_c)->tablebc

full_join(tableab, tablebc)%>%
  filter(!is.na(p_value))%>%
  pivot_wider(names_from = "group_b", values_from = "p_value")%>%
  
 
  arrange(group_a)%>%

  rename(group_name = group_a)%>%
  mutate(num = c(1:6))%>%

  relocate(group_name, num)%>%


  rename_at(vars(3:last_col()), ~as.character(c(1:6)))->results
  
  


###########################################################
#pull metadata for results object. Save in "metadata" df. 

names(pairwise)%>%
  as_tibble()%>%
  rename(col1 = value)%>%
  full_join(tibble(col1 = c("date.of.analysis")))%>%
  filter(col1 != "p.value")%>%
  arrange(col1)%>%
  mutate(col2 = c((pairwise)[[2]],
                  as.character(Sys.time()),
                  (pairwise)[[1]],
                    (pairwise)[[4]]
                  ))->metadata




##############################################################
#output

#creates google sheet with both data tables as separate sheets.
#uses Sys.Date to name sheet with today's date. 

gs4_create(paste(as.character(Sys.time()), 
                 
    #name of output file
                 name_of_file),
           
 
    
    
sheets = list(results = results, test_info = metadata))






       