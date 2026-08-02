# 2. CIBERSORTx fractions and GSVA scores
# CIBERSORTx is an external input; this file validates its table and combines
# it with locally calculated GSVA scores.
suppressPackageStartupMessages({library(GSVA); library(readxl); library(tidyverse)})

root <- normalizePath(".")
# External host-response tables and cohort manifests share one staged input tree.
input <- file.path(root, "data/analysis_inputs/host_response")
output <- file.path(root, "data/results/2_cibersortx_and_gsva")
dir.create(output, recursive=TRUE, showWarnings=FALSE)

# Produce sample annotations from the two preprocessed sample manifests
datasets <- c("GSE17156", "GSE198449")
sample_tables <- map(datasets, function(dataset) {
  samples <- read.delim(file.path(root, "data/analysis_inputs", dataset, "samples.tsv"), check.names=FALSE)
  samples$sample_id <- paste(dataset, samples$sample_id, sep="__")
  samples$dataset <- dataset
  samples$cohort <- if (dataset == "GSE198449") "SARS_CoV_2" else {
    description <- tolower(samples$local_description)
    ifelse(grepl("influenza", description), "IFV", ifelse(grepl("rhinovirus", description), "HRV", "RSV"))
  }
  samples[, c("sample_id", "dataset", "cohort", "clinical_group")]
})
metadata <- sample_tables %>% bind_rows() %>% column_to_rownames("sample_id")
# Keep this compact manifest as the shared sample-level input for step 3.
write.table(metadata, file.path(output, "sample_metadata.tsv"), sep="\t", quote=FALSE, col.names=NA)

# The CIBERSORTx table has one sample per row (Mixture as its sample identifier),
# 22 cell-fraction columns, and three diagnostic columns (P-value, Correlation, RMSE).
cibersort <- read.delim(file.path(input, "cibersortx_raw_output.tsv"), check.names=FALSE)
names(cibersort)[names(cibersort)=="Mixture"] <- "sample_id"
# Reorder the external result to the expression/metadata sample order and use
# sample IDs as row names so that the metadata join is an explicit keyed join.
cibersort <- cibersort[match(rownames(metadata), cibersort$sample_id), ]
rownames(cibersort) <- cibersort$sample_id
fraction_columns <- setdiff(names(cibersort), c("sample_id", "P-value", "Correlation", "RMSE", "Absolute score"))
# Convert fractions to numeric and check whether each sample's estimated cell composition sums approximately to one
cibersort <- cibersort %>%
  mutate(across(all_of(fraction_columns), as.numeric),
         fraction_sum=rowSums(pick(all_of(fraction_columns))),
         fraction_sum_ok=abs(fraction_sum - 1) <= 0.02)
# Append dataset and phenotype labels, producing the analysis-ready CIBERSORTx table.
cibersort <- cibersort %>% left_join(metadata %>% rownames_to_column("sample_id"), by="sample_id") %>% column_to_rownames("sample_id")
write.table(cibersort, file.path(output, "cibersortx_fractions.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# Use the already integrated expression matrix; rows are genes and columns are samples.
expression <- as.matrix(read.delim(file.path(input, "gsva_expression.tsv"), row.names=1, check.names=FALSE))
gene_set_names <- c("IFNA2", "IFNB1", "IFNW1", "IFNG", "IFN Core", "TNF")
gene_set_table <- read_excel(file.path(root, "data/analysis_inputs/reference/Supplementary Table.xlsx"), sheet="S2_Interferon-related Transcrip", skip=1)
gene_sets <- map(gene_set_table[gene_set_names], ~unique(trimws(as.character(.x[!is.na(.x)]))))
names(gene_sets) <- gene_set_names
# Keep only genes present in the expression matrix so every score uses valid features.
gene_sets <- map(gene_sets, intersect, rownames(expression))
# gsvaParam defines the score model; gsva summarizes each fixed gene set for every sample.
scores <- as.data.frame(t(as.matrix(gsva(gsvaParam(expression, gene_sets), verbose=FALSE)))) %>%
  rownames_to_column("sample_id") %>%
  left_join(metadata %>% rownames_to_column("sample_id"), by="sample_id") %>%
  column_to_rownames("sample_id")
write.table(scores, file.path(output, "gsva_scores_with_metadata.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

# Compare each feature between the two target-condition groups and the control group
# within every cohort, preserving sample counts and correcting all tests by BH.
compare_groups <- function(data, features, file) {
  comparisons <- crossing(feature=features, cohort=unique(data$cohort), group=c("asymptomatic", "symptomatic"))
  result <- pmap_dfr(comparisons, function(feature, cohort, group) {
    subset <- data %>% filter(.data$cohort == cohort)
    control <- subset %>% filter(clinical_group == "control") %>% pull(all_of(feature))
    case <- subset %>% filter(clinical_group == group) %>% pull(all_of(feature))
    tibble(feature=feature, cohort=cohort, comparison=paste0(group, "_vs_control"),
           n_control=length(control), n_group=length(case),
           pvalue=if (length(control) && length(case)) wilcox.test(case, control)$p.value else NA_real_)
  }) %>% mutate(padj=p.adjust(pvalue, method="BH"))
  write.table(result, file, sep="\t", row.names=FALSE, quote=FALSE)
}
compare_groups(cibersort, fraction_columns, file.path(output, "cibersortx_statistics.tsv"))
compare_groups(scores, names(gene_sets), file.path(output, "gsva_statistics.tsv"))
