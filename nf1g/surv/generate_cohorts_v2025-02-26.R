source("libs.R")
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(grid)
library(ggkm)
library(scales)
library(survival)
library(survminer)
library(RColorBrewer)


########dataset read in and construction######################
#reading in current dataset. 
coh1<-readRDS("nf1g/surv/cohorts-2025-01-07.rds")%>%
  mutate(aod = as.numeric(aod))

df0<-coh1%>%

  mutate(event = if_else(is.na(exp_end_date),1,0))%>%
  
  #making a convenient variable to easily compare nf1 vs. non-nf1s. 
  mutate(nf1=substr(resultant_geno, 1, 6))%>%
  mutate(geno_bg=substr(resultant_geno, 8, nchar(resultant_geno)))%>%

  mutate(aod=if_else(event==0, 150,aod))%>%
  
  #add hist category for not included in hist. categorization
  mutate(hist=if_else(is.na(hist)&event==0, "xh-ne", hist))%>%
  mutate(hist=if_else(is.na(hist)&event==1, "xh-e", hist))
  
df1<-df0%>%
  
  
  #temp selection to vis. df
  select("mouse_num", "resultant_geno", "hist", "event")%>%
  mutate(hist=toupper(hist))%>%
  
  #get 1-letter/number shorthand for histological category.
  mutate(hist_cat=substr(hist, 1,1))%>%


  #get hist_grade from hist information column. note that excluded from histology will
  #still show a "grade" but it's not actually a grade. 
  mutate(hist_grade=substr(hist, 3,4))%>%
  mutate(hist_grade=gsub("\\.", "", hist_grade))%>%
  mutate(hist_grade=gsub("15", "1.5", hist_grade))%>%
  mutate(hist_catgrade=paste0(hist_cat, ".", hist_grade))%>%
  
  
  mutate(across(c(hist_cat, hist_grade, hist_catgrade), 
               ~ if_else(hist == "NED" | grepl("^X", hist), hist, .x)))%>%
  mutate(across(c(hist_cat, hist_grade, hist_catgrade), 
                ~ if_else(hist == "NED" | grepl("^X", hist), hist, .x)))

    




nf1_renaming_dataset1 <- readxl::read_excel("nf1g/surv/nf1-renaming dataset.xlsx")
nf1_renaming_dataset2 <- readxl::read_excel("nf1g/surv/nf1-renaming dataset.xlsx", sheet = "resultant_geno")
nf1_renaming_dataset3 <- readxl::read_excel("nf1g/surv/nf1-renaming dataset.xlsx", sheet = "patho_cat")
#nf1_renaming_dataset4 <- readxl::read_excel("nf1g/surv/nf1-renaming dataset.xlsx", sheet = "patho_cat2")
nf1_renaming_dataset5 <- readxl::read_excel("nf1g/surv/nf1-renaming dataset.xlsx", sheet = "patho_cat_det")%>%
  rename(hist=patho_cat_det)

df2<-df1%>%
  left_join(nf1_renaming_dataset1)%>%
  left_join(nf1_renaming_dataset2)%>%
  left_join(nf1_renaming_dataset3)%>%
  #left_join(nf1_renaming_dataset4)
  left_join(nf1_renaming_dataset5)
  mutate(full_cohort="1")%>%
  mutate(patho_cat_name=if_else(hist=="3.1", "Pretumorigenic", patho_cat_name))



mutate(hist=if_else(!is.na(exclude_from_hist_reason)&event==1, "excluded from hist (event)", hist))%>%
  mutate(hist=if_else(!is.na(exclude_from_hist_reason)&event==0, "excluded from hist (no event)", hist))%>%



name<-df1%>%select(patho_cat_name)%>%unique()


name_order<-c("Gliomas", "Nerve sheath tumors","Spindle and epithelioid tumors",
"Glioneuronal tumors",    "Pretumorigenic",  "No histological classification", "No evidence of disease")       




df1<-df1%>%
  mutate(patho_cat_name=factor(patho_cat_name, levels=name_order))



aspectratio<-.6















genotypes<-df1%>%select(resultant_geno)%>%unique()%>%arrange()%>%pull()


geno_subset<-df1%>%
 #filter(resultant_geno==genotypes[1]|resultant_geno==genotypes[2])%>%
  select(resultant_geno, patho_cat_name)%>%
  group_by(resultant_geno,patho_cat_name)%>%
  summarise(n=n(),
            )%>%
  ungroup()%>%
  group_by(resultant_geno)%>%
  summarize(patho_cat_name, n=n,total_n=sum(n))%>%
  mutate(perc=n/total_n*100)%>%
  ungroup()%>%
  complete(resultant_geno, patho_cat_name, fill = list(perc = -1))



geno_subset


prop_test_results <- geno_subset %>%
  group_by(patho_cat_name) %>%
  summarise(
    p_value = prop.test(
      x = n,  # The count of successes (e.g., number of "green" individuals)
      n = total_n  # The total sample size
    )$p.value  # Extract p.value directly from the result
  )






