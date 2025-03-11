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
library(purrr)

source("ggplot_draw_square.R")


########dataset read in and construction######################
#reading in current dataset. 
coh1<-readRDS("nf1g/surv/cohorts-2025-03-07-v2.rds")
col_map<-readRDS("nf1g/surv/colors_map_surv.rds")



df1<-coh1%>%
#most conservative exclusion criteria
filter(is.na(exclude))%>%
  filter(!is.na(include_in_surv))%>%
  #filter(is.na(exclude_hist))%>%
  filter(resultant_geno=="nf1 KO; pten KO; ink KO; atrx KO"|resultant_geno=="nf1 KO; pten KO; ink KO; atrx wt")%>%
  mutate(resultant_geno=factor(resultant_geno, levels=c("nf1 KO; pten KO; ink KO; atrx KO","nf1 KO; pten KO; ink KO; atrx wt")))

aspectratio<-.6



# genotypes<-df1%>%select(resultant_geno)%>%unique()%>%arrange()%>%pull()



df_props<-df1%>%

  select(resultant_geno, hist_cat_name)%>%
  group_by(resultant_geno,hist_cat_name)%>%
  
  #count occurrences of each hist_cat per geno.
  reframe(n=n() )%>%
  
  #calculate percentage/likelihood of hist_cat per geno.
  group_by(resultant_geno)%>%
  reframe(hist_cat_name, n=n,total_n=sum(n))%>%
  mutate(perc=n/total_n*100)%>%

  #add in a dummy number for each hist_cat-geno combo to ensure consistency in plot display groups.
  complete(resultant_geno, hist_cat_name, fill = list(perc = -.5))%>%
  
  #get back info for resultant_geno_name ("proper" cohort names), to use if desired.
  left_join(coh1%>%select(resultant_geno, resultant_geno_name)%>%unique())%>%
  mutate(resultant_geno = factor(resultant_geno, levels = names(col_map$resultant_geno)))%>%
  mutate(resultant_geno_name = factor(resultant_geno_name, levels = names(col_map$resultant_geno_name)))%>%
  mutate(hist_cat_name = factor(hist_cat_name, levels = names(col_map$hist_cat_name)))



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
  select(hist_cat_name, p_value, estimate_WT, estimate_KO)  # Adjust columns as needed










a<-fisher.test(table(df1$resultant_geno, df1$hist_cat_name))

tumor_counts <- table(df1$resultant_geno, df1$hist_cat_name)  # Contingency table
total_counts <- rowSums(tumor_counts)  # Total tumors per genotype

# Compute proportions for each tumor type within each genotype
proportions <- sweep(tumor_counts, 1, total_counts, "/")  # Equivalent to tumor_counts / total_counts

# Compute standard error (SE) for proportions
SE <- sqrt(proportions * (1 - proportions) / total_counts)

# Convert to tidy format
prop_df <- as.data.frame(as.table(proportions)) %>%
  rename(resultant_geno = Var1, hist_cat_name = Var2, Proportion = Freq)

SE_df <- as.data.frame(as.table(SE)) %>%
  rename(resultant_geno = Var1, hist_cat_name = Var2, SE = Freq)

# Merge proportions and standard errors into one tidy tibble
result_df <- left_join(prop_df, SE_df, by = c("resultant_geno", "hist_cat_name"))




df_props1<-df_props%>%
  left_join(result_df)%>%
  janitor::clean_names()





df_props1$hist_cat_name_numeric <- as.numeric(as.factor(df_props1$hist_cat_name))*0.6
#cat_levels <- unique(df_props$hist_cat_name)
#df_props$hist_cat_name_numeric <- match(df_props$hist_cat_name, cat_levels) * 0.6  # Adjust 0.8 to fine-tune spacing



p2<-ggplot(df_props1) +
  geom_col(aes(x = hist_cat_name_numeric, 
               fill = resultant_geno,
               y = perc),
           width = 0.4,  
           position = position_dodge(0.5),
           key_glyph = draw_square)  +
  
  geom_errorbar(aes(x=hist_cat_name_numeric, 
                    ymin=perc-se*100, 
                    ymax=perc+se*100,
                    group = resultant_geno),
                position = position_dodge(0.5),
                width=.15,
                linewidth=.01,
                alpha=.5)+
  
 

  theme_classic()+

  # ggnewscale::new_scale(new_aes = "black")+
  
  
  
  scale_fill_manual(values=col_map$resultant_geno,
                    
                    
                    labels = scales::label_wrap(20))+ #wrap long strings into multiple lines.
  
  theme(
    axis.text.x = element_text(size=12,angle = 45, hjust = 1),
    plot.margin = margin(5, 5, 25, 50),
    plot.caption = element_text(hjust = 0, size = 10),
    legend.text = element_text(size = 8, hjust = 0, vjust=0.5), 
    #legend.key.height = unit(5, "mm"),
    legend.key.spacing.y = unit(5, 'mm'),
    aspect.ratio=.6)+
  
  
  ggtitle("Tumor type proportion by cohort")+
  labs(fill=NULL,
       x=NULL,
       y="% of Each Cohort",
       caption = "**Error bars represent standard error.")+
  
  
  
  
  
  
  scale_x_continuous(breaks = unique(df_props1$hist_cat_name_numeric),
                     labels = unique(df_props1$hist_cat_name))+

  
  guides(fill = guide_legend(byrow = TRUE))
  
  
  
  
  


p2




metadata_text <- paste0("src: ", 
                        rstudioapi::getSourceEditorContext()$path %>%
                          sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.) %>%
                          sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.), 
                        "\n", "time: ", 
                        lubridate::round_date(Sys.time(), "second"))

text_grob <- textGrob(metadata_text, x=.05, just="left", gp = gpar(fontsize = 9, col = "gray30"))

p2src<-grid.arrange(p2, text_grob, ncol = 1, heights = c(3, 0.3))

grid.draw(p2src)


ggsave("nf1g/surv/pub/tumor_types-perc-noXH-v0-atrxcompo_only-geno.pdf",
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       plot=p2src,
       limitsize = FALSE,
       
       
       height=10,
       width=15,
       scale = 1,
       dpi=600,
       
       
       
)  

##not necessary to plot proportions other way with only 2 groups.

# 
# stop("next plot")
# 
# df_props$resultant_geno_name_numeric <- as.numeric(as.factor(df_props$resultant_geno_name))
# geno_levels <- unique(df_props$resultant_geno_name)
# df_props$resultant_geno_name_numeric <- match(df_props$resultant_geno_name, geno_levels) * 0.6  # Adjust 0.8 to fine-tune spacing
# 
# 
# 
# 
# 
# 
# p3<-ggplot(df_props) +
#   geom_col(aes(x = resultant_geno_name_numeric, 
#                fill = hist_cat_name,
#                y = perc),
#            width = 0.4,  
#            position = position_dodge(0.5)) +  
#   theme_classic()+
#   
#   
#   
#   
#   scale_fill_manual(values=col_map$hist_cat_name)+
#   
#   theme(
#     axis.text.x = element_text(size=12,angle = 45, hjust = 1),
#     plot.margin = margin(100, 100, 100, 100),
#     aspect.ratio=.6
#   ) +
#   ggtitle("Tumor type prevalence by cohort")+
#   labs(fill=NULL,
#        x=NULL,
#        y="% tumor incidence")+
#   
#   
#   
#   scale_x_continuous(
#     breaks = unique(df_props$resultant_geno_name_numeric),
#     labels = str_wrap(unique(df_props$resultant_geno), width = 18)  # Wrap long labels
#   )+
#   
#   
#   geom_segment(aes(
#     x = as.numeric(resultant_geno_name_numeric) - 0.1,   # Slightly offset from bars
#     xend = as.numeric(resultant_geno_name_numeric) + 0.1, # Slightly offset from bars
#     y = -2.5,   # Adjust as needed for the vertical position of the line
#     yend = -2.5,  # Keep line horizontal at the same position
#     color = resultant_geno  # Color the line based on resultant_geno_name
#   ), linewidth = 1.5) +  # Line width
#   scale_color_manual(values = col_map$resultant_geno) +
#   
#   guides(color = "none") 
# 
# 
# 
# metadata_text <- paste0("src: ", 
#                         rstudioapi::getSourceEditorContext()$path %>%
#                           sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.) %>%
#                           sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.), 
#                         "\n", "time: ", 
#                         lubridate::round_date(Sys.time(), "second"))
# 
# text_grob <- textGrob(metadata_text, x=.05, just="left", gp = gpar(fontsize = 9, col = "gray30"))
# 
# p3src<-grid.arrange(p3, text_grob, ncol = 1, heights = c(3, 0.3))
# 
# 
# 
# 
# 
# 
# ggsave("nf1g/surv/pub/tumor_types-perc2-noXH-v0-atrx-comp_only-geno.pdf",
#        
#        title=paste0("src: ",
#                     
#                     rstudioapi::getSourceEditorContext()$path%>%
#                       sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
#                       sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
#                     
#                     " at ", 
#                     
#                     lubridate::round_date(Sys.time(), "second")
#        ),
#        
#        plot=p3src,
#        limitsize = FALSE,
#        
#        
#        height=10,
#        width=15,
#        scale = 1,
#        dpi=600,
#        
#        
#        
# )  
# 
