




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






