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
library(forcats)

source("ggplot_draw_square.R")
stderrorcalc <- function(count, total) {
  p <- count / total
  sqrt((p * (1 - p)) / total)
}

########dataset read in and construction######################
#reading in current dataset. 
coh1<-readRDS("nf1g/surv/cohorts-2025-03-07-v2.rds")
col_map<-readRDS("nf1g/surv/colors_map_surv.rds")



df1<-coh1%>%
#most conservative exclusion criteria
filter(is.na(exclude))%>%
  filter(!is.na(include_in_surv))%>%
  filter(is.na(exclude_hist))%>%
  filter(hist_cat_name!="No histological classification")

aspectratio<-.6



df_props0<-df1%>%

  select(resultant_geno, hist_cat_name)%>%
  group_by(resultant_geno,hist_cat_name)%>%
  
  #count occurrences of each hist_cat per geno.
  reframe(n=n() )%>%
  
  #calculate percentage/likelihood of hist_cat per geno.
  group_by(resultant_geno)%>%
  reframe(hist_cat_name, n=n,total_n=sum(n))%>%
  mutate(perc=n/total_n*100)%>%


  
  #get back info for resultant_geno_name ("proper" cohort names), to use if desired.
  left_join(coh1%>%select(resultant_geno, resultant_geno_name)%>%unique())%>%
  mutate(resultant_geno = factor(resultant_geno, levels = names(col_map$resultant_geno)))%>%
  mutate(resultant_geno_name = factor(resultant_geno_name, levels = names(col_map$resultant_geno_name)))%>%
  

  mutate(hist_cat_name = fct_drop(hist_cat_name))




group_ns<-df_props0%>%
  group_by(hist_cat_name)%>%
  reframe(n, resultant_geno, 
          hist_total_n=sum(n))

hist_cat_ns<-group_ns%>%select(hist_cat_name, hist_total_n)%>%unique()



df_props1<-df_props0%>%
  select(resultant_geno, hist_cat_name, n, perc)%>%
  filter(!is.na(n))%>%
  group_by(hist_cat_name)%>%
  reframe(norm_n=perc,
          n_actual=n,
          norm_n_total=sum(perc),
          resultant_geno)%>%
  group_by(hist_cat_name)%>%
  reframe(norm_100_n=norm_n*(100/norm_n_total), 
          hist_cat_name, n_actual, resultant_geno,
          norm_100_n_total=100)




df_props2<-df_props0%>%
  left_join(df_props1)%>%
  #add in a dummy number for each hist_cat-geno combo to ensure consistency in plot display groups.
  complete(resultant_geno, hist_cat_name, fill = list(norm_100_n = -.5))%>%
  left_join(group_ns)%>%
  mutate(x_label_ns=if_else(!is.na(n), paste0(hist_cat_name, "\n (n = ", hist_total_n, ")"), 
                            NA_character_))




df_props3<-df_props2%>%
  select(resultant_geno, hist_cat_name,
         n, total_n, 
         perc, 
         norm_100_n, 
         hist_total_n)%>%
  
  rename(total_n_geno=total_n)%>%
  rename(total_n_hist=hist_total_n)%>%
  
  select(resultant_geno,
         hist_cat_name,
         n, total_n_hist, total_n_geno)%>%
  mutate(n=if_else(is.na(n), 0, n))%>%
  
  rename(n_actual=n)%>%
  group_by(resultant_geno)%>%
  
  reframe(resultant_geno,
          hist_cat_name, n_actual,
          total_n_geno=sum(n_actual))%>%
  
  group_by(hist_cat_name)%>%
  reframe(resultant_geno,
          hist_cat_name, n_actual,
          total_n_geno,
          total_n_hist=sum(n_actual))%>%
  
  
  filter(total_n_geno!=0)%>%
  filter(total_n_hist!=0)%>%

  
  mutate(n_norm_coh_size=n_actual/total_n_geno*100)%>%
  
  group_by(hist_cat_name)%>%
  
  reframe(across(everything()), total_n_hist_norm=sum(n_norm_coh_size))%>%
  
  mutate(norm_prop=n_norm_coh_size/total_n_hist_norm)%>%
  
  mutate(norm_count=norm_prop*total_n_hist)%>%
  
  mutate(se = stderrorcalc(norm_count, total_n_hist))%>%
  
  mutate(perc=norm_prop*100)%>%
  mutate(perc=norm_count)%>%
  
  mutate(se_scaled=se*100/total_n_hist)












df_props3$hist_cat_name_numeric <- as.numeric(as.factor(df_props3$hist_cat_name))
cat_levels <- unique(df_props3$hist_cat_name)
df_props3$hist_cat_name_numeric <- match(df_props3$hist_cat_name, cat_levels) * 0.6  # Adjust 0.8 to fine-tune spacing



p2<-ggplot(df_props3) +
  geom_col(aes(x = hist_cat_name_numeric, 
               fill = resultant_geno,
               y = perc),
           width = 0.4,  
           position = position_dodge(0.5),
           key_glyph = draw_square)  +
  
  
  geom_errorbar(aes(
    x = hist_cat_name_numeric,
    ymin = ifelse(perc > 0, perc - se_scaled, NA),
    ymax = ifelse(perc > 0, perc + se_scaled, NA),
    group = resultant_geno
  ),
  position = position_dodge(0.5),
  width = .2,
  linewidth = 1,
  alpha=.5)+
  
  
  

  theme_classic()

p2


  scale_fill_manual(values=col_map$resultant_geno,
                    
                    
                   labels = scales::label_wrap(20))+ #wrap long strings into multiple lines.
  
  theme(
    axis.text.x = element_text(size=12,angle = 45, hjust = 1),
    plot.margin = margin(10, 10, 10, 20),
    plot.caption = element_text(hjust = 0, size = 10))+
  
    
    
    ggtitle("Cohort prevalence by tumor type")+
  
    
    
    labs(fill=NULL,
       x=NULL,
       y="% of Each Cohort",
       caption = "**Error bars represent standard error.")+
  
  
  
  scale_x_continuous(breaks = unique(df_props2$hist_cat_name_numeric),
                     labels = na.omit(unique(df_props2$x_label_ns)))+
  theme(
    legend.text = element_text(size = 8, hjust = 0, vjust=0.5), 
    #legend.key.height = unit(5, "mm"),
    legend.key.spacing.y = unit(5, 'mm'))+
  
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







ggsave("nf1g/surv/pub/pub_plots/tumor_types-prop-tt_cohort.pdf",
       
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


