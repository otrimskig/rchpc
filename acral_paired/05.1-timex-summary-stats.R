source("libs.R")
library(tidyverse)
library(dtplyr)
library(purrr)

#part I: generate summary stats for pathways. 


#load relevant created from gsva analysis. 
gsva_u<-readRDS("acral_paired/ds/gsva_u.rds")

standard_error <- function(x) {
  sd(x) / sqrt(length(x))
}

summary_stats_a <- gsva_u%>%as.data.frame()%>%
  rownames_to_column("pathway")%>%
  rowwise()

columns_to_select_for_stats<-colnames(summary_stats_a)[2:length(colnames(summary_stats_a))]


summary_stats<-summary_stats_a%>%
  mutate(min_value = min(c_across(all_of(columns_to_select_for_stats))),
         max_value = max(c_across(all_of(columns_to_select_for_stats))),
         mean_value = mean(c_across(all_of(columns_to_select_for_stats))),
         sd_value = sd(c_across(all_of(columns_to_select_for_stats))),
         se_value = standard_error(c_across(all_of(columns_to_select_for_stats)))
         )%>%
  select(-all_of(columns_to_select_for_stats))%>%
  ungroup()





#part II: assign samples to groups based on variable of choice. 
#Generate comparison stats for each pathway.



name_map<-readRDS("acral_paired/ds/name_map.rds")



pathways_diff<-gsva_u%>%as.data.frame()%>%
  rownames_to_column("pathway")%>%
  pivot_longer(all_of(columns_to_select_for_stats), names_to="mouse_num_type_clean", values_to="gsea_val")%>%
  left_join(name_map)%>%
  group_by(sample_type, pathway)%>%
  janitor::clean_names()



####generate p values

######################## unpaired ############

# Extract unique pathway names
pathway_names <- unique(pathways_diff$pathway)

# Initialize a tibble to store results
all_ps_unpaired <- tibble(pn = character(), pv = numeric())

# Loop through each pathway
for (pn in pathway_names) {
  # Filter data for the current pathway
  pathway_data <- pathways_diff %>%
    filter(pathway == pn)
  
  pv <- t.test(
    pathway_data %>% filter(sample_type == "subq") %>% arrange(mouse_num) %>% pull(gsea_val),
    pathway_data %>% filter(sample_type == "acral") %>% arrange(mouse_num) %>% pull(gsea_val),
    paired = FALSE
  )$p.value
  
  # Append results to the tibble
  all_ps_unpaired <- bind_rows(all_ps_unpaired, tibble(pn = pn, pv = pv))
}




######################## paired ################
########################        ################

# Extract unique pathway names
pathway_names <- unique(pathways_diff$pathway)

# Initialize a tibble to store results
all_ps_paired <- tibble(pn = character(), pv = numeric())

# Loop through each pathway
for (pn in pathway_names) {
  # Filter data for the current pathway
  pathway_data <- pathways_diff %>%
    filter(pathway == pn)
  
  pv <- t.test(
    pathway_data %>% filter(sample_type == "subq") %>% arrange(mouse_num) %>% pull(gsea_val),
    pathway_data %>% filter(sample_type == "acral") %>% arrange(mouse_num) %>% pull(gsea_val),
    paired = TRUE
  )$p.value

  # Append results to the tibble
  all_ps_paired <- bind_rows(all_ps_paired, tibble(pn = pn, pv = pv))
}


##################################################

##################################################
#generate p-value corrections. 
unpaired_corrected <- all_ps_unpaired %>%
  mutate(
    pv_bf = p.adjust(pv, method = "bonferroni"),
    pv_fdr = p.adjust(pv, method = "fdr")
  ) %>%
  rename_with(~ paste0(., "_up"), .cols = pv:last_col())


paired_corrected <- all_ps_paired %>%
  mutate(
    pv_bf = p.adjust(pv, method = "bonferroni"),
    pv_fdr = p.adjust(pv, method = "fdr")
  )%>%
  rename_with(~ paste0(., "_pa"), .cols = pv:last_col())


corrected_ps<-left_join(unpaired_corrected, paired_corrected)%>%
  mutate(diff=pv_up-pv_pa)




############################################
#generate all summary stats

all_stats<-pathways_diff%>%
  summarise(mean=mean(gsea_val),
            sd=sd(gsea_val),
            n=n())%>%
  ungroup()%>%
  pivot_wider(names_from = "sample_type", values_from = c("mean", "sd", "n"))%>%
  janitor::clean_names()%>%
  
  mutate(fc=mean_acral/mean_subq)%>%
  
  left_join(corrected_ps%>%rename(pathway=pn))%>%
  select(-diff)%>%
  group_by(pathway)%>%
  mutate(n=n())%>%
  
  left_join(summary_stats)




saveRDS(all_stats, "acral_paired/ds/gsva_pathway_stats.rds")



