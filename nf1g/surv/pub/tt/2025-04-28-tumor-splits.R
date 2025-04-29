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


########dataset read in and construction######################
#reading in current dataset. 
coh1<-readRDS("nf1g/surv/cohorts-2025-03-07-v2.rds")
col_map<-readRDS("nf1g/surv/colors_map_surv.rds")


df1<-coh1%>%
  #most conservative exclusion criteria
  filter(is.na(exclude))%>%
  filter(!is.na(include_in_surv))


prop_var_name<-"hist_grade_name"


df_props<-df1%>%
  
  select(resultant_geno, !!sym(prop_var_name))%>%
  group_by(resultant_geno, !!sym(prop_var_name))%>%
  
  #count occurrences of each hist_cat per geno.
  reframe(n=n() )%>%
  
  #calculate percentage/likelihood of hist_cat per geno.
  group_by(resultant_geno)%>%
  reframe(!!sym(prop_var_name), n=n,total_n=sum(n))%>%
  mutate(perc=n/total_n*100)%>%
  
  #add in a dummy number for each hist_cat-geno combo to ensure consistency in plot display groups.
  complete(resultant_geno,!!sym(prop_var_name), fill = list(perc = -.5))%>%
  
  #get back info for resultant_geno_name ("proper" cohort names), to use if desired.
  left_join(coh1%>%select(resultant_geno, resultant_geno_name)%>%unique())%>%
  mutate(resultant_geno = factor(resultant_geno, levels = names(col_map$resultant_geno)))%>%
  mutate(resultant_geno_name = factor(resultant_geno_name, levels = names(col_map$resultant_geno_name)))




df_props_w <- df_props %>%
  select(-perc, -resultant_geno_name)%>%
  pivot_wider(
    names_from = resultant_geno,    # Pivot by `resultant_geno`
    values_from = c(n, total_n),    # Get `n` and `total_n` columns
    names_glue = "{resultant_geno}_{.value}"  # Create custom column names like 'genotype_n' and 'genotype_total_n'
  )



results <- df_props_w %>%
  rowwise() %>%  # Apply row-wise operations
  mutate(
    prop_test = list({
      # Check if any count is less than a threshold (e.g., 5)
      if(any(c(`nf1 KO; pten KO; ink KO; atrx KO_n`, `nf1 KO; pten KO; ink KO; atrx wt_n`) < 5)) {
        # Use Fisher's Exact Test for small counts
        fisher.test(matrix(c(`nf1 KO; pten KO; ink KO; atrx KO_n`, `nf1 KO; pten KO; ink KO; atrx wt_n`, 
                             `nf1 KO; pten KO; ink KO; atrx KO_total_n`, `nf1 KO; pten KO; ink KO; atrx wt_total_n`), 
                           nrow = 2))
      } else {
        # Use prop.test() otherwise
        prop.test(
          x = c(`nf1 KO; pten KO; ink KO; atrx KO_n`, `nf1 KO; pten KO; ink KO; atrx wt_n`),  # Success counts (x)
          n = c(`nf1 KO; pten KO; ink KO; atrx KO_total_n`, `nf1 KO; pten KO; ink KO; atrx wt_total_n`)  # Total counts (n)
        )
      }
    }),
    p_value = prop_test$p.value,
    estimate_WT = prop_test$estimate[1],  # Estimate for WT genotype
    estimate_KO = prop_test$estimate[2]   # Estimate for KO genotype
  ) %>%
  select(!!sym(prop_var_name), p_value, estimate_WT, estimate_KO)  # Adjust columns as needed




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

# Merge proportions and standard errors into one tidy tibble
result_df <- left_join(prop_df, SE_df, by = c("resultant_geno", setNames(prop_var_name, prop_var_name)))





#######plot 5 prep#############

df_props_pen0 <- df_props %>%
  mutate(tumor = if_else(hist_grade_name %in% c("No evidence of disease", "Excluded from histology (no event)"),
                         "No tumor", "Tumor"))%>%
  # 
  mutate(tumor = if_else(hist_grade_name %in% c("Pretumorigenic"),
                         "Pretumorigenic", tumor))%>%

  mutate(perc=if_else(perc<0, NA, perc))%>%
  filter(!is.na(perc))%>%
  
  mutate(resultant_geno_label=str_wrap(resultant_geno, width=8))%>%
  
  mutate(tumor=factor(tumor, levels=c("Tumor", "Pretumorigenic", "No tumor")))

#######plot5 start########################################################################### 


p5<-ggplot(df_props_pen0, aes(x = tumor, y = perc, fill=hist_grade_name)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9, color="black")+
  facet_grid(~resultant_geno_label)+
  
 scale_fill_manual(values=col_map$hist_grade_name)+
  
  
  theme_pubclean()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        guides(fill = guide_legend(nrow = 9)))+  # Adjust size as needed
  labs(x = NULL,
       y="% of cohort")







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

