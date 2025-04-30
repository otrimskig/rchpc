source("libs.R")
suppressMessages(source("nf1g/colors_map_create.R", echo = FALSE))

library(tidyverse)
library(ggpubr)
library(gridExtra)
library(ggkm)
library(scales)
library(survival)
library(survminer)
library(RColorBrewer)

source("ggplot_draw_square.R")


######## dataset read in and construction######################
#reading in current dataset. 
coh1<-readRDS("nf1g/surv/cohorts-2025-03-07-v2.rds")
col_map<-readRDS("nf1g/surv/colors_map_surv.rds")


df1<-coh1%>%
  #most conservative exclusion criteria
  filter(is.na(exclude))%>%
  filter(!is.na(include_in_surv))


prop_var_name<-"hist_grade_name"


df_props<-df1%>% #prep dataset and count proportions of each prop_var_name per genotype.
  
  select(resultant_geno, !!sym(prop_var_name))%>%
  group_by(resultant_geno, !!sym(prop_var_name))%>%
  
  #count occurrences of each hist_cat per geno.
  reframe(n=n() )%>%
  
  #calculate percentage/likelihood of hist_cat per geno.
  group_by(resultant_geno)%>%
  reframe(!!sym(prop_var_name), n=n)%>%
  
  #complete categories. add n=0 for any non-occurring combinations.
  complete(resultant_geno,!!sym(prop_var_name), fill = list(n = 0))%>%
  
  group_by(resultant_geno)%>%
  reframe(!!sym(prop_var_name), n=n, total_n=sum(n))%>%
  
  mutate(perc=n/total_n*100)%>%
  #get back info for resultant_geno_name ("proper" cohort names), to use if desired.
  left_join(coh1%>%select(resultant_geno, resultant_geno_name)%>%unique())%>%
  mutate(resultant_geno = factor(resultant_geno, levels = names(col_map$resultant_geno)))%>%
  mutate(resultant_geno_name = factor(resultant_geno_name, levels = names(col_map$resultant_geno_name)))%>%
  
  mutate(across(where(is.factor), ~ droplevels(.)))

  
  #add in a dummy number for each hist_cat-geno combo to ensure consistency in plot display groups.
  #should be done after stats, before plotting.
  #complete(resultant_geno,!!sym(prop_var_name), fill = list(perc = -.5))%>%
  
 


# 
# 
# df_props_w <- df_props %>%
#   select(-perc, -resultant_geno_name)%>%
#   pivot_wider(
#     names_from = resultant_geno,    # Pivot by `resultant_geno`
#     values_from = c(n, total_n),    # Get `n` and `total_n` columns
#     names_glue = "{resultant_geno}_{.value}"  # Create custom column names like 'genotype_n' and 'genotype_total_n'
#   )
# 


# results <- df_props_w %>%
#   rowwise() %>%  # Apply row-wise operations
#   mutate(
#     prop_test = list({
#       # Check if any count is less than a threshold (e.g., 5)
#       if(any(c(`nf1 KO; pten KO; ink KO; atrx KO_n`, `nf1 KO; pten KO; ink KO; atrx wt_n`) < 5)) {
#         # Use Fisher's Exact Test for small counts
#         fisher.test(matrix(c(`nf1 KO; pten KO; ink KO; atrx KO_n`, `nf1 KO; pten KO; ink KO; atrx wt_n`, 
#                              `nf1 KO; pten KO; ink KO; atrx KO_total_n`, `nf1 KO; pten KO; ink KO; atrx wt_total_n`), 
#                            nrow = 2))
#       } else {
#         # Use prop.test() otherwise
#         prop.test(
#           x = c(`nf1 KO; pten KO; ink KO; atrx KO_n`, `nf1 KO; pten KO; ink KO; atrx wt_n`),  # Success counts (x)
#           n = c(`nf1 KO; pten KO; ink KO; atrx KO_total_n`, `nf1 KO; pten KO; ink KO; atrx wt_total_n`)  # Total counts (n)
#         )
#       }
#     }),
#     p_value = prop_test$p.value,
#     estimate_WT = prop_test$estimate[1],  # Estimate for WT genotype
#     estimate_KO = prop_test$estimate[2]   # Estimate for KO genotype
#   ) %>%
#   select(!!sym(prop_var_name), p_value, estimate_WT, estimate_KO)  # Adjust columns as needed




tumor_counts <- table(df1$resultant_geno,  df1[[prop_var_name]])  # Contingency table
total_counts <- rowSums(tumor_counts)  # Total tumors per genotype

# Compute proportions for each tumor type within each genotype
proportions <- sweep(tumor_counts, 1, total_counts, "/")  # Equivalent to tumor_counts / total_counts



# Compute standard error (SE) for proportions
SE <- sqrt(proportions * (1 - proportions) / total_counts)

# Convert to tidy format
prop_df <- as.data.frame(as.table(proportions)) %>%
  rename(resultant_geno = Var1,   !!prop_var_name := Var2, Proportion = Freq)

SE_df <- as.data.frame(as.table(SE)) %>%
  rename(resultant_geno = Var1,  !!prop_var_name := Var2, SE = Freq)

# Merge proportions and standard errors into one tidy tibble back into original proportion dataset.
df_props1 <- left_join(prop_df, SE_df, by = c("resultant_geno", setNames(prop_var_name, prop_var_name)))%>%
  rename(proportion=Proportion, se_prop=SE)%>%
  left_join(df_props)



stats_list<-list(meta=list(src=rstudioapi::getSourceEditorContext()$path%>%
                             sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                             sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                           
                           datetime=lubridate::round_date(Sys.time(), "second")),
                 
                 full_props=df_props1)




####### 3g prep#############

df_props_3g <- df_props %>%
  mutate(tumor = if_else(hist_grade_name %in% c("No evidence of disease", "Excluded from histology (no event)"),
                         "No tumor", "Tumor"))%>%
  # 
  mutate(tumor = if_else(hist_grade_name %in% c("Pretumorigenic"),
                         "Pretumorigenic", tumor))%>%

  mutate(perc=if_else(perc<0, NA, perc))%>%
  filter(!is.na(perc))%>%
  
  mutate(resultant_geno_label=str_wrap(resultant_geno, width=8))%>%
  
  mutate(tumor=factor(tumor, levels=c("Tumor", "Pretumorigenic", "No tumor")))


############## 3g stats #############

prop_tests_df<-df_props_3g%>%
  select(resultant_geno, hist_grade_name, n, tumor)%>%
  group_by(resultant_geno, tumor)%>%
  summarise(ncount=sum(n))%>%
  pivot_wider(values_from = ncount,  names_from = tumor)%>%
  janitor::clean_names()%>%
  mutate(total_n=tumor+no_tumor+pretumorigenic)



# All pairwise 2x3 comparisons
pairwise_tests3g <- combn(1:nrow(prop_tests_df), 2, simplify = FALSE) %>%
  map_df(~{
    g1 <- prop_tests_df[.x[1], ]
    g2 <- prop_tests_df[.x[2], ]
    
    mat <- rbind(
      c(g1$tumor, g1$pretumorigenic, g1$no_tumor),
      c(g2$tumor, g2$pretumorigenic, g2$no_tumor)
    )
    
    fisher_p <- fisher.test(mat, simulate.p.value = TRUE, B = 1e5)$p.value
    chi_p <- chisq.test(mat, correct = FALSE)$p.value
    
    tibble(
      group1 = as.character(g1$resultant_geno),
      group2 = as.character(g2$resultant_geno),
      fisher_p = fisher_p,
      chisq_p = chi_p
    )
  }) %>%
  arrange(group1, group2)





####### 2g prep#############

df_props_2g <- df_props %>%
  mutate(tumor = if_else(hist_grade_name %in% c("No evidence of disease", "Pretumorigenic", "Excluded from histology (no event)"),
                         "No tumor", "Tumor"))%>%

  mutate(perc=if_else(perc<0, NA, perc))%>%
  filter(!is.na(perc))%>%
  
  mutate(resultant_geno_label=str_wrap(resultant_geno, width=8))%>%
  
  mutate(tumor=factor(tumor, levels=c("Tumor", "No tumor")))




############## 3g stats #############

prop_tests_df<-df_props_2g%>%
  select(resultant_geno, hist_grade_name, n, tumor)%>%
  group_by(resultant_geno, tumor)%>%
  summarise(ncount=sum(n))%>%
  pivot_wider(values_from = ncount,  names_from = tumor)%>%
  janitor::clean_names()%>%
  mutate(total_n=tumor+no_tumor)



# All pairwise 2x3 comparisons
pairwise_tests2g <- combn(1:nrow(prop_tests_df), 2, simplify = FALSE) %>%
  map_df(~{
    g1 <- prop_tests_df[.x[1], ]
    g2 <- prop_tests_df[.x[2], ]
    
    mat <- rbind(
      c(g1$tumor,g1$no_tumor),
      c(g2$tumor, g2$no_tumor)
    )
    
    fisher_p <- fisher.test(mat, simulate.p.value = TRUE, B = 1e5)$p.value
    chi_p <- chisq.test(mat, correct = FALSE)$p.value
    
    tibble(
      group1 = as.character(g1$resultant_geno),
      group2 = as.character(g2$resultant_geno),
      fisher_p = fisher_p,
      chisq_p = chi_p
    )
  }) %>%
  arrange(group1, group2)





#######################################












table_data <- with(df_props1, tapply(n, list(hist_grade_name, resultant_geno), sum))
chisq.test(table_data)



library(nnet)
df_props1$hist_grade_name <- relevel(df_props1$hist_grade_name, ref = "No evidence of disease")
model <- multinom(hist_grade_name ~ resultant_geno, weights = n, data = df_props1)
summary(model)
















####### plot5 start########################################################################### 


p5<-ggplot(df_props_pen0, aes(x = tumor, y = perc, fill=hist_grade_name)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9, color="black")+
  facet_grid(~resultant_geno_label)+
  
 scale_fill_manual(values=col_map$hist_grade_name)+
  
  
  theme_pubclean()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        guides(fill = guide_legend(nrow = 9)))+  # Adjust size as needed
  labs(x = NULL,
       y="% of cohort")






####### save plot5 start####
metadata_text <- paste0("src: ", 
                        rstudioapi::getSourceEditorContext()$path %>%
                          sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.) %>%
                          sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.), 
                        "\n", "time: ", 
                        lubridate::round_date(Sys.time(), "second"))

text_grob <- textGrob(metadata_text, x=.05, just="left", gp = gpar(fontsize = 9, col = "gray30"))

p5src<-grid.arrange(p5, text_grob, ncol = 1, heights = c(3, 0.3))




ggsave("nf1g/surv/pub/pub_plots/penetrance-stacked-grade0.pdf",
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       plot=p5src,
       limitsize = FALSE,
       
       
       height=3,
       width=5,
       scale = 4,
       dpi=600,
       
       
       
)  
ggsave("nf1g/surv/pub/pub_plots/penetrance-stacked-grade0b.pdf",
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       plot=p5src,
       limitsize = FALSE,
       
       
       height=5,
       width=3.5,
       scale = 2,
       dpi=600,
       
       
       
)  


