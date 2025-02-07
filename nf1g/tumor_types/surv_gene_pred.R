source("libs.R")
library(tidyverse)
library(survival)
library(caret)


per_sample<-readRDS("nf1g/ds/v10-per_sample_updated.rds")%>%
  select(sample_id, aod)



cohorts2<-readRDS("nf1g/surv/cohorts-2025-01-07.rds")%>%
  filter(!is.na(rna_seq_sample_id))%>%
  dplyr::select(aod, event, mouse_num, rna_seq_sample_id)%>%
  rename(sample_id=rna_seq_sample_id)




all_info<-readRDS("nf1g/ds/vm-02-filtered_rpkms.rds")%>%
  select(gene_name_ms, sample_id, rpkm)%>%
  left_join(cohorts2)

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

dim(expression_mat)


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
  p_value = sapply(cox_results, function(res) {
    if (!is.null(res)) res$coefficients[1, "Pr(>|z|)"] else NA
  }),
  hazard_ratio = sapply(cox_results, function(res) {
    if (!is.null(res)) exp(res$coefficients[1, "coef"]) else NA
  })
)


saveRDS(results_df, "nf1g/tumor_types/aod-gene-analysis.rds")



