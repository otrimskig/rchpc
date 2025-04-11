source("libs.R")
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(gridtext)
library(grid)
library(gridExtra)

library(circlize)
library(purrr)
library(tidyverse)
library(ggplot2)
library(crayon)
library(fs)



mat0<-readRDS("acral_sub_rppa/ds/umat-per_sample_375_rel_diff_th.rds")
stats0<-readRDS("acral_sub_rppa/ds/stats_fr_threshold_var.rds")
sample0<-readRDS("acral_sub_rppa/ds/sample_info0.rds")




sig_abs<-stats0%>%
  filter(t_test_p_excl<=0.05)%>%
  pull(antibody_name)



mat1<-mat0[sig_abs,]





hm1<-Heatmap(zmat1)
g_hm1<-grid.grabExpr(draw(draw((hm1), annotation_legend_side = "left")))


ggsave("acral_sub_rppa/plots/hm-filtered01.pdf",
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       plot=g_hm1,
       limitsize = FALSE,
       
       
       height=6,
       width=3,
       scale = 3
)



umat1<-mat1

zmat1<-t(scale(t(mat1)))




zdf0 <- zmat1 %>%
  as_tibble(rownames = "antibody_name") %>%
  pivot_longer(cols = -antibody_name,
               names_to = "sample_num_type",
               values_to = "rppa_val") %>%
  separate(sample_num_type, into = c("sample_num", "sample_type"), sep = "_")


z_stats <- zdf0 %>%
  group_by(antibody_name, sample_type) %>%
  summarise(
    mean = mean(rppa_val, na.rm = TRUE),
    sd = sd(rppa_val, na.rm = TRUE),
    se = sd / sqrt(sum(!is.na(rppa_val))),  # Corrected to exclude NA
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = sample_type,
    values_from = c(mean, sd, se),
    names_glue = "{.value}_{sample_type}"
  )

zdf1<-zdf0%>%
  left_join(z_stats)%>%
  mutate(diff=mean_FP-mean_SQ)%>%
  mutate(updown=if_else(diff>0,  "up", "down"))%>%
  mutate(updown=factor(updown, levels=c("up", "down")))




# Plot
ggplot(data = zdf1, aes(
  x = rppa_val,
  y = fct_reorder(antibody_name, diff),  # Corrected y-axis ordering
  color = sample_type
)) +
  geom_segment(aes(
    x = mean_FP-se_FP, xend = mean_FP+se_FP,
    y = antibody_name, yend = antibody_name
  ),color ="#F8766D", alpha = 0.02, linewidth = 1.5) +
  
  
  geom_segment(aes(
    x = mean_SQ-se_SQ, xend = mean_SQ+se_SQ,
    y = antibody_name, yend = antibody_name
   
  ),color="#00BFC4",alpha = 0.02, linewidth = 1.5) +
  
  
  
 #geom_point(aes(x = mean_FP, y = antibody_name), shape = 124, size = 4, color = "#F8766D") +
  #geom_point(aes(x = mean_SQ, y = antibody_name), shape = 124, size = 4, color = "#00BFC4") +
  
  geom_point(size = 3, alpha = 0.5, position = position_dodge(width=2)) +
  
  theme_bw() +
  theme(
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "lines"),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    legend.justification = "top",
    legend.direction = "vertical"
  )



library(ggplot2)
library(forcats)
library(ggstance)  # for vertical dodge

p_dot1<-ggplot(data = zdf1, aes(
  x = rppa_val,
  y = fct_reorder(antibody_name, diff),
  color = sample_type
)) +
  # geom_segment(aes(
  #   x = mean_FP - se_FP, xend = mean_FP + se_FP,
  #   y = antibody_name, yend = antibody_name
  # ),
  # color = "#F8766D", alpha = 0.02, linewidth = 1.5
  # 
  # ) +
  # geom_segment(aes(
  #   x = mean_SQ - se_SQ, xend = mean_SQ + se_SQ,
  #   y = antibody_name, yend = antibody_name
  # ),
  # color = "#00BFC4", alpha = 0.02, linewidth = 1.5,
  # position = position_dodgev(height = 0.5)
  # ) +
  # 
  
facet_grid(vars(updown), scales = "free_y", space="free")+
  
  
  geom_point(size = 5, alpha = 0.7, shape=1)+
  geom_point(size = 5, alpha = .3,shape=16)+
  
  
  theme_bw() +
  theme(
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "lines"),
    axis.text.y = element_text(size = 13),
    legend.position = "right",
    legend.justification = "top",
    legend.direction = "vertical"
  )



p_dot1









ggsave("acral_sub_rppa/plots/dotplot-filtered01.pdf",
       
       title=paste0("src: ",
                    
                    rstudioapi::getSourceEditorContext()$path%>%
                      sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                      sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                    
                    " at ", 
                    
                    lubridate::round_date(Sys.time(), "second")
       ),
       
       plot=p_dot1,
       limitsize = FALSE,
       
       
       height=6,
       width=3,
       scale = 2.5
)












udf0 <- umat1 %>%
  as_tibble(rownames = "antibody_name") %>%
  pivot_longer(cols = -antibody_name,
               names_to = "sample_num_type",
               values_to = "rppa_val") %>%
  separate(sample_num_type, into = c("sample_num", "sample_type"), sep = "_")



unscaled_means <- udf0 %>%
  group_by(antibody_name, sample_type) %>%
  summarise(mean = mean(rppa_val, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = sample_type, values_from = mean) %>%
  rename_with(.cols = -antibody_name, .fn = ~paste0("mean_", .))


udf1<-udf0%>%
  left_join(unscaled_means)%>%
  mutate(diff=mean_FP-mean_SQ)




# Plot
ggplot(data = udf1, aes(
  x = rppa_val,
  y = fct_reorder(antibody_name, diff),  # Corrected y-axis ordering
  color = sample_type
)) +
  geom_segment(aes(
    x = mean_FP, xend = mean_SQ,
    y = antibody_name, yend = antibody_name
  ),
  color = "gray", alpha = 0.1, linewidth = 2) +
  
  geom_point(aes(x = mean_FP, y = antibody_name), shape = 124, size = 4, color = "#E41A1C") +
  geom_point(aes(x = mean_SQ, y = antibody_name), shape = 124, size = 4, color = "#377EB8") +
  
  geom_point(size = 3, alpha = 0.35) +
  
  theme_bw() +
  theme(
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "lines"),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    legend.justification = "top",
    legend.direction = "vertical"
  )




