
# Merge proportions and standard errors into one tidy tibble back into original proportion dataset.
df_props1 <- left_join(prop_df, SE_df, by = c("resultant_geno", setNames(prop_var_name, prop_var_name)))%>%
  rename(proportion=Proportion, se_prop=SE)%>%
  left_join(df_props)%>%
  
  mutate(resultant_geno_label=str_wrap(resultant_geno,
                                       width = 17))

ggplot(df_props1, aes(x = hist_grade_name, y = perc, fill=hist_grade_name)) +
  geom_bar(stat = "identity", width = 0.9, color="black") +
  facet_grid(~resultant_geno_label) +
  scale_fill_manual(values=col_map$hist_grade_name) +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        guides(fill = guide_legend(nrow = 9))) +
  labs(x = NULL, y="% of cohort") +
  geom_text(aes(label = round(perc, 1)), vjust = -0.5, size = 4)  # Adjust vjust for positioning and size as needed




ggplot(df_props1, aes(x = hist_grade_name, y = perc, fill = hist_grade_name)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  facet_grid(~resultant_geno_label) +
  scale_fill_manual(values = col_map$hist_grade_name) +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        guides(fill = guide_legend(nrow = 9))) +
  labs(x = NULL, y = "% of cohort") +
  geom_text(aes(label = ifelse(perc > 0, round(perc, 1), "")), 
            vjust = -0.5, size = 4, 
            position = position_nudge(y = 0.1))  # Adjust vertical nudging to avoid overlap


ggplot(df_props1, aes(x = hist_grade_name, y = perc, fill = hist_grade_name)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  facet_grid(~resultant_geno_label) +
  scale_fill_manual(values = col_map$hist_grade_name) +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        guides(fill = guide_legend(nrow = 9))) +
  labs(x = NULL, y = "% of cohort") +
  
  
  geom_label(aes(label = ifelse(perc > 0, round(perc, 1), "")), 
             vjust = -0.5, size = 4, 
             position = position_nudge(y = 0.1),
             label.padding = unit(0.2, "lines"),  # Padding around the text
             label.background = "white",  # White background
             label.border = element_rect(color = "black", size = 0.5))  # Black border








df_props1a<-df_props1%>%
  mutate(label_text=if_else(perc!=0, paste0(sprintf("%.1f", perc), " (",  n,  ")"), NA))%>%
  mutate(resultant_geno_label=paste0(resultant_geno_label, "\n (", total_n, ")"))




ggplot(df_props1a, aes(x = hist_grade_name, y = perc, fill = hist_grade_name)) +
  geom_bar(stat = "identity", width = 0.9, color = "black") +
  facet_grid(~resultant_geno_label) +
  scale_fill_manual(values = col_map$hist_grade_name) +
 
  #theme_pubclean() +
  #theme_bw()+
  #theme_classic2()+
  
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        guides(fill = guide_legend(nrow = 9))) +
  labs(x = NULL, y = "% of cohort") +
  geom_text(aes(label = label_text), 
            #vjust = -0.5, 
            hjust=-.1,
            size = 4, 
            position = position_nudge(y = 0.1)) +
  
  ylim(c(0,100))+
  
  labs(hist_grade_name="")+
  
  coord_flip()

