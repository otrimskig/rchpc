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
library(ggh4x) 

source("ggplot_draw_square.R")


add_ws <- function(label, n_final_string) {
  n_ws <- max(n_final_string - nchar(label), 0) %/% 2  # Ensure non-negative and even distribution
  paste0(strrep(" ", n_ws), label)
}



########dataset read in and construction######################
#reading in current dataset. 
coh1<-readRDS("nf1g/surv/cohorts-2025-03-07-v2.rds")
col_map<-readRDS("nf1g/surv/colors_map_surv.rds")
aspectratio<-.6

df1<-coh1%>%
  #most conservative exclusion criteria
  filter(is.na(exclude))%>%
  filter(!is.na(include_in_surv))%>%
  filter(hist_cat_name!="No histological classification")
  #
df1%>%
  select(starts_with("hist"))%>%
  filter(hist=="PT"|hist=="3.1")%>%
  distinct()


  
#filter(is.na(exclude_hist))
  #filter(as.numeric(hist_grade)>1)%>%
  

ct_geno<-df1%>%
  group_by(resultant_geno)%>%
  
  summarise(count=n())


ct_tt<-df1%>%
  group_by(hist_cat_name)%>%
  
  summarise(count=n())

sum_df1<-df1%>%
  group_by(resultant_geno, hist_cat_name, hist_grade_name)%>%
  
  summarise(n_geno=n())%>%
  
  ungroup()%>%
  
  
  
  arrange(resultant_geno, hist_cat_name, hist_grade_name)
 
  
   
  


comp_1<-sum_df1%>%
  filter(grepl("\\d$", hist_grade_name))%>%
  mutate(across(where(is.factor), droplevels))%>%
  complete(resultant_geno, hist_cat_name, hist_grade_name,
           fill = list(n_geno = 0))


comp_2<-sum_df1%>%
  filter(!grepl("\\d$", hist_grade_name))%>%
  #filter(!grepl("No grade assigned", hist_grade_name))
  mutate(across(where(is.factor), droplevels))%>%

  complete(resultant_geno, hist_cat_name,
           fill = list(n_geno = 0))%>%
  
  mutate(hist_grade_name=if_else(hist_cat_name=="Glioneuronal tumors"|
                                   hist_cat_name=="No histological classification", "No grade assigned", hist_grade_name))%>%
  
  mutate(hist_grade_name=if_else(is.na(hist_grade_name), hist_cat_name, hist_grade_name))






comp_3<-full_join(comp_1, comp_2)%>%
  mutate(hist_grade_name = factor(hist_grade_name, levels = names(col_map$hist_grade_name)))












sum_df2<-comp_3%>%
  ungroup%>%group_by(resultant_geno)%>%
  reframe(across(everything()),
          n_total_geno=sum(n_geno),
          prop_geno=n_geno/n_total_geno,
          se_prop = sqrt((prop_geno * (1 - prop_geno)) / n_total_geno))%>%
  
  group_by(hist_cat_name)%>%
  
          
  reframe(across(everything()), n_total_hist=sum(n_geno),
            prop_total_hist=sum(prop_geno),
          prop_prop_total_hist=prop_geno/prop_total_hist,
          se_prop_scaled=se_prop/prop_total_hist)%>%

  relocate(resultant_geno)%>%
  relocate(n_total_hist, .after=n_total_geno)
  
  
  






# cat_levels <- unique(sum_df2$hist_grade_name)
# sum_df2$hist_grade_name_numeric <- match(sum_df2$hist_grade_name_numeric, cat_levels) * 0.2  # Adjust 0.8 to fine-tune spacing




#plot data part 1

gg_data1<-sum_df2%>%
  filter(grepl("\\d$", hist_grade_name)|
           hist_grade_name=="No grade assigned")%>%
  #          hist_grade_name=="No evidence of disease")%>%
  #filter(hist_cat_name!="No histological classification")%>%
  
  mutate(across(where(is.factor), droplevels))%>%
  
  complete(hist_grade_name, hist_cat_name, resultant_geno, fill=list(n_geno=0))
  
  
 
  

gg_data2<-sum_df2%>%
  filter(!grepl("\\d$", hist_grade_name))%>%
  
  filter(hist_grade_name=="No evidence of disease")%>%
  filter(hist_cat_name!="No histological classification"&hist_grade_name!="No grade assigned")%>%
  
  mutate(across(where(is.factor), droplevels))%>%
  
  complete(hist_grade_name, hist_cat_name, resultant_geno, fill=list(n_geno=0))



gg_data3<-sum_df2%>%
  filter(!grepl("\\d$", hist_grade_name))%>%
  
  filter(hist_grade_name=="Pretumorigenic")%>%
  filter(hist_cat_name!="No histological classification"&hist_grade_name!="No grade assigned")%>%
  
  mutate(across(where(is.factor), droplevels))%>%
  
  complete(hist_grade_name, hist_cat_name, resultant_geno, fill=list(n_geno=0))
  
gg_data4<-sum_df2%>%
  filter(grepl("^Excluded", hist_cat_name))%>%
  
  mutate(across(where(is.factor), droplevels))%>%
  
  complete(hist_grade_name, resultant_geno, fill=list(n_geno=0))





gg_data5<-bind_rows(gg_data2, gg_data3, gg_data4)%>%
  mutate(hist_grade_name="No tumor or uncharacterized")%>%
  mutate(across(where(is.factor), droplevels))
  








gg_data<-bind_rows(gg_data1, gg_data5)%>%
  filter(grepl("\\d$",hist_grade_name))%>%
  mutate(across(where(is.factor), droplevels))%>%
  
  mutate(hist_grade_name_numeric=as.numeric(factor(hist_grade_name))*.3)%>%
  mutate(prop_geno_dummy=if_else(n_geno==0, -.5, prop_geno*100))






#plot part 1
tt_facet<-
ggplot(gg_data) +
  geom_col(aes(x = hist_cat_name,
               fill = resultant_geno,
               y = prop_geno_dummy,
               alpha = ifelse(prop_geno_dummy < 0, .4, NA)
              ),
           #color = if_else(prop_geno_dummy < 0, "black", "transparent"),
           width = .8,
           position = position_dodge(width = .85,  preserve = "single"),
           key_glyph = draw_square)+
  
  scale_fill_manual(values=col_map$resultant_geno)+
  scale_alpha_identity(guide = "none")+
  theme_classic()+
  
  facet_wrap(vars(hist_grade_name), ncol=1)+
  
  geom_errorbar(aes(
    x = hist_cat_name,
    ymin = ifelse(prop_geno_dummy > 0, prop_geno_dummy-se_prop * 100, NA),
    ymax = ifelse(prop_geno_dummy > 0, prop_geno_dummy+se_prop * 100, NA),
    group=resultant_geno
  ),
  position = position_dodge(width = .85,  preserve = "single"),
  width = 0.3, linewidth = 0.2,
  alpha=.5)+
  
 
  ggtitle("Penetrance by grade and tumor type")+
  
  labs(fill=NULL,
       x=NULL,
       y="% of Each Cohort",
       caption = "**Error bars represent standard error of proportion.")+
  # 

  theme(
  axis.text.x = element_text(size=12,angle = 45, hjust = 1),
  strip.text = element_text(size = 10),  # Increase the size of facet labels
  plot.margin = margin(10, 50, 10, 50),
  plot.caption = element_text(hjust = 0, size = 10),
  axis.line.x.bottom = element_line(color = "black", size = 0.5),
  axis.line.y.left = element_line(color = "black", size = 0.5),
  panel.grid.major = element_blank(),  # Removes major grid lines
  panel.grid.minor = element_blank(),
  legend.position = "none" )





# tt_facet$labels$title


#save plot
plot_object_name_for_session<-"tt_facet"
plot_filename_output<-"nf1g/surv/pub/pub_plots/tt_facet_3.pdf"


metadata_text <- paste0("src: ", 
                        rstudioapi::getSourceEditorContext()$path %>%
                          sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.) %>%
                          sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.), 
                        "\n", "time: ", 
                        lubridate::round_date(Sys.time(), "second"))

text_grob <- grid::textGrob(metadata_text, x=.05, just="left", gp = grid::gpar(fontsize = 9, col = "gray30"))

p3src<-gridExtra::grid.arrange(get(plot_object_name_for_session), text_grob, ncol = 1, heights = c(3, 0.3))

ggsave(filename=plot_filename_output,
       
       plot=p3src,
       limitsize = FALSE,
       
       height=9,
       width=6,
       scale = 1,
       dpi=600,
       
       title=paste0(get(plot_object_name_for_session)$labels$title,
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ))
       
       

