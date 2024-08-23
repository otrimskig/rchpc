source("libs.R")
library(tidyverse)
library(dtplyr)
library(purrr)


#load relevant created from gsva analysis. 
gsva_u<-readRDS("k19mf/ds/gsva_u-onco.rds")

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
  ungroup()%>%
  column_to_rownames("pathway")%>%
  as.matrix.data.frame()

saveRDS(summary_stats, "k19mf/ds/gsva_pathway_stats-onco.rds")



grouping<-readRDS("k19mf/ds/vm-00-sample_info.rds")%>%
  select(mouse_num, tumor_type)


pathways_diff<-gsva_u%>%as.data.frame()%>%
  rownames_to_column("pathway")%>%
  pivot_longer(all_of(columns_to_select_for_stats), names_to="mouse_num", values_to="gsea_val")%>%
  left_join(grouping)%>%
  group_by(tumor_type, pathway)%>%
  janitor::clean_names()
  
pathways_diff2<-pathways_diff%>%
  summarise(mean=mean(gsea_val),
          sd=sd(gsea_val),
          n=n())%>%
  ungroup()%>%
  pivot_wider(names_from = "tumor_type", values_from = c("mean", "sd", "n"))%>%
  janitor::clean_names()%>%

  mutate(fc=mean_foot_pad/mean_subq)




pathway_names<-unique(pathways_diff$pathway)
all_ps<-tibble(pn=NA, pv=NA)

for (i in 1:length(pathway_names)){


pn<-pathway_names[i]

x<-pathways_diff%>%filter(pathway==pathway_names[[i]])%>%filter(tumor_type=="subq")%>%pull(gsea_val)

y<-pathways_diff%>%filter(pathway==pathway_names[[i]])%>%filter(tumor_type=="foot pad")%>%pull(gsea_val)

pv<-t.test(x, y,
       alternative = "two.sided",
       paired = FALSE)$p.value


all_ps<-bind_rows(all_ps, tibble(pn, pv))

}

all_ps<-all_ps%>%
  filter(!is.na(pv))
