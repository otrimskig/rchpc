source("libs.R")
library(tidyverse)
library(survival)
library(caret)



per_sample_all<-readRDS("nf1g/ds/v10-per_sample_updated.rds")

tumor_types<-per_sample_all%>%
  select(patho_cat_name)%>%
  unique()%>%
  pull()%>%
  c()%>%
  as.character()

tumor_types<-tumor_types[1:3]






per_sample<-readRDS("nf1g/ds/v10-per_sample_updated.rds")%>%
  select(sample_id, aod, patho_cat_name)


cohort_info<-readRDS("nf1g/surv/cohorts-2025-01-07.rds")%>%
  filter(!is.na(rna_seq_sample_id))%>%
  dplyr::select(aod, event, mouse_num, rna_seq_sample_id)%>%
  rename(sample_id=rna_seq_sample_id)%>%
  mutate(aod=as.numeric(aod))%>%
  
  left_join(per_sample)%>%
  filter(!is.na(patho_cat_name))






for (i in 1:length(tumor_types)){
  
spec_samples<-cohort_info%>%filter(patho_cat_name==tumor_types[i])  
  


all_info<-readRDS("nf1g/ds/vm-02-filtered_rpkms.rds")%>%
  select(gene_name_ms, sample_id, rpkm)%>%
  semi_join(spec_samples)%>%
  left_join(spec_samples)

sample_surv<-all_info%>%
  select(sample_id, aod, event)%>%
  unique()%>%
  mutate(aod=as.numeric(aod))%>%
  mutate(event=as.numeric(event))



expression_mat<-all_info%>%
  select(gene_name_ms, sample_id, rpkm)%>%
  pivot_wider(names_from = sample_id, values_from = rpkm)%>%
  column_to_rownames("gene_name_ms")%>%
  as.matrix.data.frame()




 # For nearZeroVar
nzv_genes <- nearZeroVar(expression_mat)

expression_mat <- scale(expression_mat)

cox_results <- apply(expression_mat, 1, function(gene) {
  tryCatch({
    coxph(Surv(aod, event) ~ gene, data = sample_surv) %>% summary()
  }, error = function(e) {
    message("Skipping gene due to error: ", e$message)
    return(NULL)  # Skip the gene
  })
})


results_df <- data.frame(
  gene = rownames(expression_mat),
  coef = sapply(cox_results, function(res) {
    if (!is.null(res)) res$coefficients[1, "coef"] else NA
  }),
  hazard_ratio = sapply(cox_results, function(res) {
    if (!is.null(res)) exp(res$coefficients[1, "coef"]) else NA
  }),
  se_coef = sapply(cox_results, function(res) {
    if (!is.null(res)) res$coefficients[1, "se(coef)"] else NA
  }),
  z_value = sapply(cox_results, function(res) {
    if (!is.null(res)) res$coefficients[1, "z"] else NA
  }),
  p_value = sapply(cox_results, function(res) {
    if (!is.null(res)) res$coefficients[1, "Pr(>|z|)"] else NA
  }),
  lower_95_CI = sapply(cox_results, function(res) {
    if (!is.null(res)) res$conf.int[1, "lower .95"] else NA
  }),
  upper_95_CI = sapply(cox_results, function(res) {
    if (!is.null(res)) res$conf.int[1, "upper .95"] else NA
  }),
  concordance = sapply(cox_results, function(res) {
    if (!is.null(res)) res$concordance[1] else NA
  }),
  concordance_se = sapply(cox_results, function(res) {
    if (!is.null(res)) res$concordance[2] else NA
  }),
  likelihood_ratio_p = sapply(cox_results, function(res) {
    if (!is.null(res)) res$logtest["pvalue"] else NA
  }),
  wald_test_p = sapply(cox_results, function(res) {
    if (!is.null(res)) res$waldtest["pvalue"] else NA
  }),
  score_logrank_p = sapply(cox_results, function(res) {
    if (!is.null(res)) res$sctest["pvalue"] else NA
  })
)


saveRDS(results_df, 
        paste0("nf1g/tumor_types/aod-gene-analysis-",
               tumor_types[i],
               ".rds"))



}



column_explanations <- data.frame(
  Column_Name = c("coef", "hazard_ratio", "se_coef", "z_value", "p_value", 
                  "lower_95_CI", "upper_95_CI", "concordance", "concordance_se", 
                  "likelihood_ratio_p", "wald_test_p", "score_logrank_p"),
  Description = c(
    "Estimated regression coefficient (log hazard ratio)",
    "Exponentiated coefficient (Hazard Ratio, HR)",
    "Standard error of the coefficient",
    "Wald test z-score",
    "p-value from the Wald test",
    "Lower bound of 95% confidence interval for HR",
    "Upper bound of 95% confidence interval for HR",
    "Concordance index (C-index), measures model predictive ability",
    "Standard error of the C-index",
    "p-value for the likelihood ratio test",
    "p-value for the Wald test",
    "p-value for the Score (log-rank) test"
  )
)

write_csv(column_explanations, "nf1g/tumor_types/gene_surv_analysis_stats_explanation.csv")
