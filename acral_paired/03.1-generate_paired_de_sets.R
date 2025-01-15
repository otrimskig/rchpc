##########################
source("libs.R")
library(tidyverse)
library(edgeR)
library(dtplyr)

sample_info <- readRDS("acral_paired/ds/v00-sample_info.rds")
read_counts <- readRDS("acral_paired/ds/v02-filtered_rpkms.rds")

all_data <- read_counts %>%
left_join(sample_info) %>%
ungroup() %>%
arrange(sample_id)

# View all available categories.
colnames(all_data)

# Use to determine category for selection of comparison.
categories <- c("sample_type")

r_folder_name <- "acral_paired"

# Generate DE sets from all categories specified above.
# tryCatch({
# for (i in 1:length(categories)) {

i<-1

category <- categories[i]

comp_el <- all_data %>%
ungroup() %>%
select(!!sym(category)) %>%
count(!!sym(category)) %>%
select(1) %>%
pull()

comb_mat <- combn(comp_el, 2)

for (c in 1:ncol(comb_mat)) {
ga <- comb_mat[1, c]
gb <- comb_mat[2, c]

ct <- sym(category)

comp_info <- all_data %>%
group_by(sample_id) %>% slice(1) %>% ungroup() %>%
select(sample_id:last_col()) %>%
filter(!!ct == ga | !!ct == gb)

group_a_count <- comp_info %>% filter(!!ct == ga) %>% count() %>% pull()
group_b_count <- comp_info %>% filter(!!ct == gb) %>% count() %>% pull()
}







comp_counts <- all_data %>% 
  semi_join(comp_info, by = "sample_id") %>%
  select(gene_name_ms, gene_id_ms, read_count, sample_id) %>%
  arrange(gene_name_ms)%>%
  pivot_wider(names_from = "sample_id", values_from = "read_count") %>%
  filter(!is.na(gene_name_ms)) %>%
  group_by(gene_name_ms) %>% slice(1) %>% ungroup() %>%
  group_by(gene_id_ms) %>% slice(1) %>% ungroup() %>%
  column_to_rownames("gene_name_ms")


c_counts <- comp_counts %>% select(-gene_id_ms)
c_genes <- comp_counts %>% select(gene_id_ms)

c_group <- comp_info %>%
  select(all_of(ct)) %>%
  mutate(model_bin = if_else(!!ct == ga, 0, 1)) %>%
  pull(model_bin)

c_pairs <- comp_info %>%
  select(mouse_num) %>% pull()

dge <- DGEList(counts = c_counts, genes = c_genes, group = c_group)

keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)

design <- model.matrix(~ c_pairs + c_group)

dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef = 2)

top <- topTags(lrt, n = 100000)

output1 <- as.data.frame(top) %>%
  rownames_to_column("gene_name_ms") %>%
  as_tibble()

output2 <- output1 %>%
  left_join(all_data %>%
              select(gene_name_ms, rpkm, sample_id) %>% 
              semi_join(comp_info, by = "sample_id")
            ) %>%
  pivot_wider(values_from = rpkm, names_from = sample_id, names_prefix = "rpkm_")


dir.create(paste0(r_folder_name, "/dexps"), showWarnings = FALSE)

saveRDS(output2, paste0(r_folder_name, "/dexps/", "dexp-paired", fs::path_sanitize(paste0(category, "-", ga, " v. ", gb)), ".rds"))
 



