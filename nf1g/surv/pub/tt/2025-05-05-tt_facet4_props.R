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
  mutate(prop_geno_dummy=if_else(n_geno==0, -.5, prop_geno*100))%>%
  mutate(hist_grade_num=as.character(substr(hist_grade_name, nchar(hist_grade_name), nchar(hist_grade_name))))%>%
  mutate(hist_cat_name=factor(hist_cat_name, levels=c("Gliomas", "Spindle and epithelioid tumors", "Nerve sheath tumors", "Glioneuronal tumors")))%>%
  
  mutate(hist_grade_name_dummy=if_else(prop_geno_dummy<0, paste0(hist_grade_name, "_no fill"), hist_grade_name))%>%
  
  
  
  mutate(hist_cat_name_dummy=if_else(prop_geno_dummy<0, paste0(hist_cat_name, "_no fill"), hist_cat_name))%>%
  mutate(cat_grade_combo=paste(hist_cat_name, hist_grade_name))%>%
  
  mutate(cat_grade_combo_dummy=if_else(prop_geno_dummy<0, paste0(cat_grade_combo, "_no fill"), cat_grade_combo))

  
col_list$hist_grade_name

col_list$hist_cat_name

combos<-gg_data%>%
  select(hist_cat_name, hist_grade_name, cat_grade_combo)%>%
  mutate(cat_grade_combo=paste(hist_cat_name, hist_grade_name))%>%
  unique()%>%
  mutate(hex=NA)

#plot part 1
tt_facet<-
ggplot(gg_data) +
  geom_col(aes(x = hist_cat_name,
               fill = hist_grade_name_dummy,
               y = prop_geno_dummy,
               #color = if_else(prop_geno_dummy < 0, "1", "0")
              # alpha = ifelse(prop_geno_dummy < 0, .4, NA)
              ),
          # 
           width = .8,
           position = position_dodge(width = .85,  preserve = "single"),
           key_glyph = draw_square)+
  
  scale_fill_manual(values = col_list$hist_grade_name)+
  scale_color_manual(values=c("1"="grey", "0"="transparent"))+
  
  
  scale_alpha_identity(guide = "none")+
  theme_classic()+
  
 facet_wrap(vars(resultant_geno), ncol=1)+
  
  # geom_errorbar(aes(
  #   x = hist_cat_name,
  #   ymin = ifelse(prop_geno_dummy > 0, prop_geno_dummy-se_prop * 100, NA),
  #   ymax = ifelse(prop_geno_dummy > 0, prop_geno_dummy+se_prop * 100, NA),
  #   group=resultant_geno
  # ),
  # position = position_dodge(width = .85,  preserve = "single"),
  # width = 0.3, linewidth = 0.2,
  # alpha=.5)+
  # 
 
  ggtitle("Penetrance by grade and tumor type")+
  ylim(c(-1,25))+
  labs(fill=NULL,
       x=NULL,
       y="% of Each Cohort",
       caption = "**Error bars represent standard error of proportion.")+


  theme(
  axis.text.x = element_text(size=12,angle = 45, hjust = 1),
  strip.text = element_text(size = 10),  # Increase the size of facet labels
  plot.margin = margin(10, 50, 10, 50),
  plot.caption = element_text(hjust = 0, size = 10),
  axis.line.x.bottom = element_line(color = "black", size = 0.5),
  axis.line.y.left = element_line(color = "black", size = 0.5),
  panel.grid.major = element_blank(),  # Removes major grid lines
  panel.grid.minor = element_blank())
 # legend.position = "none" )


# 














tt_facet








# #save plot
# plot_object_name_for_session<-"tt_facet"
# plot_filename_output<-"nf1g/surv/pub/pub_plots/tt_facet_3.pdf"
# 
# 
# metadata_text <- paste0("src: ", 
#                         rstudioapi::getSourceEditorContext()$path %>%
#                           sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.) %>%
#                           sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.), 
#                         "\n", "time: ", 
#                         lubridate::round_date(Sys.time(), "second"))
# 
# text_grob <- grid::textGrob(metadata_text, x=.05, just="left", gp = grid::gpar(fontsize = 9, col = "gray30"))
# 
# p3src<-gridExtra::grid.arrange(get(plot_object_name_for_session), text_grob, ncol = 1, heights = c(3, 0.3))
# 
# ggsave(filename=plot_filename_output,
#        
#        plot=p3src,
#        limitsize = FALSE,
#        
#        height=9,
#        width=6,
#        scale = 1,
#        dpi=600,
#        
#        title=paste0(get(plot_object_name_for_session)$labels$title,
#                     
#                     " at ", 
#                     
#                     lubridate::round_date(Sys.time(), "second")
#        ))
#        










#plot part 1
tt_facet2<-
  ggplot(gg_data) +
  geom_col(aes(x = hist_grade_name,
               fill = hist_cat_name,
               y = prop_geno_dummy,
               #color = if_else(prop_geno_dummy < 0, "1", "0")
               alpha = ifelse(prop_geno_dummy < 0, .4, NA)
  ),
  # 
  width = .8,
  position = position_dodge(width = .85,  preserve = "single"),
  key_glyph = draw_square)+
  
  scale_fill_manual(values = col_list$hist_cat_name)+
  scale_color_manual(values=c("1"="grey", "0"="transparent"))+
  
  
  scale_alpha_identity(guide = "none")+
  theme_classic()+
  
  facet_wrap(vars(resultant_geno), ncol=1)+
  
  # geom_errorbar(aes(
  #   x = hist_cat_name,
  #   ymin = ifelse(prop_geno_dummy > 0, prop_geno_dummy-se_prop * 100, NA),
  #   ymax = ifelse(prop_geno_dummy > 0, prop_geno_dummy+se_prop * 100, NA),
  #   group=resultant_geno
  # ),
  # position = position_dodge(width = .85,  preserve = "single"),
  # width = 0.3, linewidth = 0.2,
  # alpha=.5)+
  # 

ggtitle("Penetrance by grade and tumor type")+
  ylim(c(-1,25))+
  labs(fill=NULL,
       x=NULL,
       y="% of Each Cohort",
       caption = "**Error bars represent standard error of proportion.")+
  
  
  theme(
    axis.text.x = element_text(size=12,angle = 45, hjust = 1),
    strip.text = element_text(size = 10),  # Increase the size of facet labels
    plot.margin = margin(10, 50, 10, 50),
    plot.caption = element_text(hjust = 0, size = 10),
    axis.line.x.bottom = element_line(color = "black", size = 0.5),
    axis.line.y.left = element_line(color = "black", size = 0.5),
    panel.grid.major = element_blank(),  # Removes major grid lines
    panel.grid.minor = element_blank())
# legend.position = "none" )


# 

tt_facet2







library(colorspace)
# Define base colors
base_colors <- c(
  "Nerve sheath tumors" = "#B74C4D",
  "Spindle and epithelioid tumors" = "#5A8BA6",
  "Gliomas" = "#66B68F",
  "Glioneuronal tumors" = "#A678A9"
)

# Function to create a palette of 4 shades (light to dark) from a base color
get_grade_shades <- function(hex) {
  lighten(hex, amount = c(.9, 0.2, -.2, -.6))
   # c(hex, darken(hex, amount = 0.8))
}

# Apply to each base color
grade_colors <- purrr::map_df(names(base_colors), function(name) {
  shades <- get_grade_shades(base_colors[[name]])
  tibble(
    hist_cat_name = name,
    hist_grade_name = paste("Grade", 1:4),
    cat_grade_combo = paste(name, paste("Grade", 1:4)),
    hex = shades
  )
})

grade_colors

named_color_list <- setNames(grade_colors$hex, grade_colors$cat_grade_combo)
named_color_list













cat_exists<-gg_data%>%
  select(cat_grade_combo, n_geno)%>%
  group_by(cat_grade_combo)%>%
  summarise(exists=max(n_geno))%>%
  filter(exists>0)%>%
  select(cat_grade_combo)
  
  
  







gg_data00<-gg_data%>%
  semi_join(cat_exists)



#plot part 1
tt_facet3<-
  ggplot(gg_data00) +
  geom_col(aes(x = hist_cat_name,
               fill = cat_grade_combo_dummy,
               y = prop_geno_dummy,
               #color = if_else(prop_geno_dummy < 0, "1", "0")
               alpha = ifelse(prop_geno_dummy < 0, .4, NA)
  ),
  width = .8,
  position = position_dodge(width = .85, preserve="single"),
  key_glyph = draw_square)+
  
  scale_fill_manual(values = named_color_list)+
  scale_color_manual(values=c("1"="grey", "0"="transparent"))+
  
  scale_alpha_identity(guide = "none")+
  theme_classic()+
  
  facet_wrap(vars(resultant_geno), ncol=1)+

ggtitle("Penetrance by grade and tumor type")+
  ylim(c(-1,25))+
  labs(fill=NULL,
       x=NULL,
       y="% of Each Cohort",
       caption = "**Error bars represent standard error of proportion.")+
  
  theme(
    axis.text.x = element_text(size=12,angle = 45, hjust = 1),
    strip.text = element_text(size = 10),  # Increase the size of facet labels
    plot.margin = margin(10, 50, 10, 50),
    plot.caption = element_text(hjust = 0, size = 10),
    axis.line.x.bottom = element_line(color = "black", size = 0.5),
    axis.line.y.left = element_line(color = "black", size = 0.5),
    panel.grid.major = element_blank(),  # Removes major grid lines
    panel.grid.minor = element_blank())





tt_facet3




library(dplyr)

gg_data01 <- gg_data00 %>%
  group_by(hist_cat_name) %>%
  mutate(
    bar_within_group = row_number() - 1,
    group_index = cur_group_id(),
    x_plot = group_index * 10 + bar_within_group
  ) %>%
  ungroup()



x_breaks <- gg_data00 %>%
  group_by(hist_cat_name) %>%
  summarize(center = mean(x_plot), .groups = "drop")






tt_facet4 <-
  ggplot(gg_data01) +
  geom_col(aes(
    x = x_plot,
    fill = cat_grade_combo_dummy,
    y = prop_geno_dummy,
    alpha = ifelse(prop_geno_dummy < 0, .4, NA)
  ),
  width = 0.8,
  key_glyph = draw_square) +
  
  scale_fill_manual(values = named_color_list) +
  scale_alpha_identity(guide = "none") +
  theme_classic() +
  
  # Optional: remove facets if you now encode everything on the x-axis
  # facet_wrap(vars(resultant_geno), ncol = 1) +
  
  ggtitle("Penetrance by grade and tumor type") +
  ylim(c(-1, 25)) +
  labs(
    fill = NULL,
    x = NULL,
    y = "% of Each Cohort",
    caption = "**Error bars represent standard error of proportion."
  ) +
  
  # Add custom x-axis labels at group centers
  scale_x_continuous(
    breaks = gg_data00 %>%
      group_by(hist_cat_name) %>%
      summarize(center = mean(x_plot)) %>%
      pull(center),
    labels = gg_data00 %>%
      group_by(hist_cat_name) %>%
      summarize(label = first(hist_cat_name)) %>%
      pull(label)
  ) +
  
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    strip.text = element_text(size = 10),
    plot.margin = margin(10, 50, 10, 50),
    plot.caption = element_text(hjust = 0, size = 10),
    axis.line.x.bottom = element_line(color = "black", size = 0.5),
    axis.line.y.left = element_line(color = "black", size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )













#plot part 1
tt_facet5<-
  ggplot(gg_data) +
  geom_col(aes(x = hist_cat_name,
               fill = cat_grade_combo_dummy,
               y = prop_geno_dummy,
               #color = if_else(prop_geno_dummy < 0, "1", "0")
               alpha = ifelse(prop_geno_dummy < 0, .4, NA)
  ),
  width = .8,
  position = position_dodge(width = .85, preserve="single"),
  key_glyph = draw_square)+
  
  scale_fill_manual(values = named_color_list)+
  scale_color_manual(values=c("1"="grey", "0"="transparent"))+
  
  scale_alpha_identity(guide = "none")+
  theme_classic()+
  
  facet_wrap(vars(resultant_geno), ncol=1)+
  
  ggtitle("Penetrance by grade and tumor type")+
  ylim(c(-1,25))+
  labs(fill=NULL,
       x=NULL,
       y="% of Each Cohort",
       caption = "**Error bars represent standard error of proportion.")+
  
  theme(
    axis.text.x = element_text(size=12,angle = 45, hjust = 1),
    strip.text = element_text(size = 10),  # Increase the size of facet labels
    plot.margin = margin(10, 50, 10, 50),
    plot.caption = element_text(hjust = 0, size = 10),
    axis.line.x.bottom = element_line(color = "black", size = 0.5),
    axis.line.y.left = element_line(color = "black", size = 0.5),
    panel.grid.major = element_blank(),  # Removes major grid lines
    panel.grid.minor = element_blank())





tt_facet5

