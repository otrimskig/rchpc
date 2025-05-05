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


######## dataset read in and construction######################
#reading in current dataset. 
coh1<-readRDS("nf1g/surv/cohorts-2025-03-07-v2.rds")
col_map<-readRDS("nf1g/surv/colors_map_surv.rds")


df1<-coh1%>%
  #most conservative exclusion criteria
  filter(is.na(exclude))%>%
  filter(!is.na(include_in_surv))


prop_var_name<-"hist_grade_name"

###initial dataset prep ##################
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



## 3group (Tumor/No tumor/pretumor) ####################
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
    
    suppressWarnings(
      chi_p <- chisq.test(mat, correct = FALSE)$p.value
      )
    
    tibble(
      group1 = as.character(g1$resultant_geno),
      group2 = as.character(g2$resultant_geno),
      fisher_p = fisher_p,
      chisq_p = chi_p
    )
  }) %>%
  arrange(group1, group2)

####### plot5 start########################################################################### 


p3g<-ggplot(df_props_3g, aes(x = tumor, y = perc, fill=hist_grade_name)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9, color="black")+
  facet_grid(~resultant_geno_label)+
  
  scale_fill_manual(values=col_map$hist_grade_name)+
  
  
  theme_pubclean()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        guides(fill = guide_legend(nrow = 9)))+  # Adjust size as needed
  labs(x = NULL,
       y="% of cohort")




labels<-df_props_3g%>%
  left_join(tibble(
    resultant_geno = names(col_map$resultant_geno),
    color_hex = col_map$resultant_geno
  )
  )%>%
  select(resultant_geno_label, color_hex)%>%
  unique()


# 
# facet_levels <- levels(df_props_3g$resultant_geno)
# facet_colors <- col_map$resultant_geno[facet_levels]


strip_elements <- Map(function(fill) {
  element_rect(fill = fill, colour = "black", linewidth = 0.5)
}, labels$color_hex)


# Build plot
ggplot(df_props_3g, aes(x = tumor, y = perc, fill = hist_grade_name)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9, color = "black") +
  
  # Apply in strip_themed
  facet_grid2(~resultant_geno_label,
              strip = strip_themed(
                background_x = strip_elements,
                text_x = element_text(size = 8)
              ))+
  
  scale_fill_manual(values = col_list$hist_grade_name) +
  
  theme_pubclean() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    strip.text.x = element_text(size = 11, face = "bold", color = "black")
  ) +
  labs(x = NULL, y = "% of cohort")






col_map$resultant_geno


df_props_3g$resultant_geno%>%unique()




setdiff(levels(df_props_3g$resultant_geno), names(col_map$resultant_geno))

####### save plot5 start####
metadata_text <- paste0("src: ", 
                        rstudioapi::getSourceEditorContext()$path %>%
                          sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.) %>%
                          sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.), 
                        "\n", "time: ", 
                        lubridate::round_date(Sys.time(), "second"))

text_grob <- textGrob(metadata_text, x=.05, just="left", gp = gpar(fontsize = 9, col = "gray30"))

p3g_gr<-grid.arrange(p3g, text_grob, ncol = 1, heights = c(3, 0.3))




ggsave("nf1g/surv/pub/pub_plots/penetrance-stacked-grade0.pdf",
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       plot=p3g_gr,
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
       
       plot=p3g_gr,
       limitsize = FALSE,
       
       
       height=5,
       width=3.5,
       scale = 2,
       dpi=600,
       
       
       
)  


###### end of plot saves #########################
## 2group (Tumor/No) ###################
####### 2g prep#############

df_props_2g <- df_props %>%
  mutate(tumor = if_else(hist_grade_name %in% c("No evidence of disease", "Pretumorigenic", "Excluded from histology (no event)"),
                         "No tumor", "Tumor"))%>%

  mutate(perc=if_else(perc<0, NA, perc))%>%
  filter(!is.na(perc))%>%
  
  mutate(resultant_geno_label=str_wrap(resultant_geno, width=8))%>%
  
  mutate(tumor=factor(tumor, levels=c("Tumor", "No tumor")))


############## 2g stats #############

prop_tests_df<-df_props_2g%>%
  select(resultant_geno, hist_grade_name, n, tumor)%>%
  group_by(resultant_geno, tumor)%>%
  summarise(ncount=sum(n))%>%
  pivot_wider(values_from = ncount,  names_from = tumor)%>%
  janitor::clean_names()%>%
  mutate(total_n=tumor+no_tumor)



# All pairwise 2x2 comparisons
pairwise_tests2g <- combn(1:nrow(prop_tests_df), 2, simplify = FALSE) %>%
  map_df(~{
    g1 <- prop_tests_df[.x[1], ]
    g2 <- prop_tests_df[.x[2], ]
    
    mat <- rbind(
      c(g1$tumor,g1$no_tumor),
      c(g2$tumor, g2$no_tumor)
    )
    
    fisher_p <- fisher.test(mat, simulate.p.value = TRUE, B = 1e5)$p.value
    chi_p <- suppressWarnings(
      chisq.test(mat, correct = FALSE)$p.value
    )
    tibble(
      group1 = as.character(g1$resultant_geno),
      group2 = as.character(g2$resultant_geno),
      fisher_p = fisher_p,
      chisq_p = chi_p
    )
  }) %>%
  arrange(group1, group2)



####### 2g plot start########################################################################### 


p2g<-ggplot(df_props_2g, aes(x = tumor, y = perc, fill=hist_grade_name)) +
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

p2g_gr<-grid.arrange(p2g, text_grob, ncol = 1, heights = c(3, 0.3))




ggsave("nf1g/surv/pub/pub_plots/penetrance-stacked-grade0-2g.pdf",
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       plot=p2g_gr,
       limitsize = FALSE,
       
       
       height=3,
       width=5,
       scale = 4,
       dpi=600,
       
       
       
)  
ggsave("nf1g/surv/pub/pub_plots/penetrance-stacked-grade0-2g-b.pdf",
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       plot=p2g_gr,
       limitsize = FALSE,
       
       
       height=5,
       width=3.5,
       scale = 2,
       dpi=600,
       
       
       
)  


###### end of plot saves #########################



df_props_stack0<-df_props_2g%>%
  mutate(hist_grade_name=factor(hist_grade_name, 
                                
                                levels=rev(c("Excluded from histology (event)", 
                                                          "No grade assigned", 
                                         "Grade 4", "Grade 3", "Grade 2", "Grade 1",
                                         "Pretumorigenic", "No evidence of disease",
                                         "Excluded from histology (no event)"))))

df_stack_pen<-df_props_stack0%>%
  group_by(resultant_geno)%>%
  filter(tumor=="Tumor")%>%
  summarise(pen=sum(perc))%>%
  left_join(col_df)




df_props_stack1<-df_props_stack0%>%
  left_join(df_stack_pen)%>%
  left_join(col_df)



pstack<-ggplot(df_props_stack1, aes(x =resultant_geno, y = perc, fill=hist_grade_name)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, #color="black", 
           alpha=.7)+
  #facet_grid(~resultant_geno_label)+
  
  scale_fill_manual(values=col_map$hist_grade_name)+
  
  # 
  # geom_tile(data=df_stack_pen, aes(x=resultant_geno,
  #                                  y=-1    #pen
  #                                  ),
  #           fill=df_stack_pen$hex_col,
  #           height=.8,
  #           width=.8,
  #           #fill="black",
  #           #color="black",
  #           alpha=1)+

  
  
  geom_tile(data=df_stack_pen, aes(x=resultant_geno, 
                                   y=pen
  ), 
  fill="black", 
  height=.2,
  width=.8,
  #fill="black", 
  #color="black", 
  alpha=1)+
  
  
  
  
  
  
 
  geom_label(
    aes(y = pen, label = paste0(round(pen), "%")),
    fill = "white",
    color = NA,                 # no border
    label.r = unit(0, "pt"),    # square corners
    size = 3
  ) +
  
  # Black text overlaid
  geom_text(
    aes(y = pen, label = paste0(round(pen), "%")),
    color = "black",
    size = 3
  )+
 theme_pubr()+
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, 
                                   size = 12),
        guides(fill = guide_legend(nrow = 9)))+  # Adjust size as needed
  labs(x = NULL,
       y="% of cohort")+
  scale_y_continuous(expand = expansion(mult = c(0.0, 0.0)))+
  coord_flip()
 


pstack
# 
# table_data <- with(df_props1, tapply(n, list(hist_grade_name, resultant_geno), sum))
# chisq.test(table_data)
# 
# 
# 
# library(nnet)
# df_props1$hist_grade_name <- relevel(df_props1$hist_grade_name, ref = "No evidence of disease")
# model <- multinom(hist_grade_name ~ resultant_geno, weights = n, data = df_props1)
# summary(model)
# 







library(ggplot2)
library(patchwork)  # for plot layout

# Ensure the order of genotypes (if needed)
df_props_stack1$resultant_geno <- factor(df_props_stack1$resultant_geno)
df_stack_pen$resultant_geno <- factor(df_stack_pen$resultant_geno, 
                                      levels = levels(df_props_stack1$resultant_geno))

pstack2 <- ggplot(df_props_stack1, aes(y = resultant_geno, 
                                       x = perc, 
                                       fill = hist_grade_name)) +
  
  # Main stacked bars, now horizontal by default
  geom_bar(stat = "identity", position = "stack", width=.7, alpha = 0.7) +
  
  scale_fill_manual(values = col_map$hist_grade_name, name = NULL)+
  
  # Tile annotations now placed at x = pen
  geom_tile(
    data = df_stack_pen,
    aes(y = resultant_geno, x = pen),
    inherit.aes = FALSE,
    fill = "black",
    height = 0.9,
    width = 0.2,
    alpha = 1
  ) +
  
  
  # White label background
  geom_label(
    data = df_stack_pen,
    aes(y = resultant_geno, x = pen+2.1, label = paste0(round(pen), "%")),
    inherit.aes = FALSE,
    fill = "white",
    color = NA,
    label.r = unit(0, "pt"),
    size = 3
  ) +
  scale_y_discrete(limits = geno_order)+
  
  
  # Black overlay text
  geom_text(
    data = df_stack_pen,
    aes(y = resultant_geno, x = pen+2.1, label = paste0(round(pen), "%")),
    inherit.aes = FALSE,
    color = "black",
    size = 3
  ) +
  
  
  
  theme_pubr() +
  theme(
    axis.text.y =  element_blank(),
    plot.margin = margin(5, 50, 5, 0),
    axis.text.x = element_text(size = 12),
    guides(fill = guide_legend(nrow = 9))
  ) +
  
  labs(x = "% of cohort", y = NULL) +
  
  scale_x_continuous(expand = expansion(mult = c(0, 0)))


pstack2


geno_order <- levels(df_props_stack1$resultant_geno)

df_anno <- data.frame(
  resultant_geno = factor(geno_order, levels = geno_order),  # Explicitly factor with correct order
  y = geno_order,
  x = 1,
  color = df_stack_pen$hex_col  # assuming it's in the same order
)


# Annotation plot (one row per genotype)
p_anno <- ggplot(df_anno, aes(y = resultant_geno)) +
  
  geom_tile(aes(x = 9.9, fill = color), width = .2, height = 0.8) +
  
  geom_text(aes(x = 9.6, label = resultant_geno), hjust = 1, size = 4) +
  
  scale_fill_identity() +
  scale_y_discrete(limits = geno_order)+
  #coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    #plot.margin = margin(5, 10, 5, 5),
    axis.text = element_blank()
  ) +
  xlim(0, 10)



plot_save_loc<-"nf1g/surv/pub/pub_plots/penetrance-stacked-grade-horiz-colcoh-1.pdf"


#save plot 7
metadata_text <- paste0("save loc: ", plot_save_loc, "\n", 
  
  
  "src: ", 
                        rstudioapi::getSourceEditorContext()$path %>%
                          sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.) %>%
                          sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.), 
                        "\n", "time: ", 
                        lubridate::round_date(Sys.time(), "second"))


combined_plot7 <- p_anno + pstack2 + 
  plot_layout(widths = c(1, 3)) +
  plot_annotation(caption = metadata_text, theme = theme(plot.caption = element_text(size = 9, hjust = 0, color = "gray30")))

# Then just print it
combined_plot7



ggsave(plot_save_loc,
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       plot=combined_plot7,
       limitsize = FALSE,
       
       
       height=2.5,
       width=6,
       scale = 2,
       dpi=600,
       
       
       
)  


