# 3. COCONUT integration, WGCNA, and DEG intersection
suppressPackageStartupMessages({library(COCONUT); library(WGCNA); library(tidyverse)})
options(stringsAsFactors=FALSE)

root <- normalizePath(".")
input <- file.path(root, "data/results/2_cibersortx_and_gsva")
raw_input <- file.path(root, "data/analysis_inputs")
deg_input <- file.path(root, "data/results/1_deg_and_enrichment")
output <- file.path(root, "data/results/3_coconut_wgcna_intersection")
dir.create(output, recursive=TRUE, showWarnings=FALSE)

# Step 2 exports the shared metadata; select the original sample20 training set.
datasets <- c("GSE17156", "GSE198449")
samples <- read.delim(file.path(input, "sample_metadata.tsv"), row.names=1, check.names=FALSE)
set.seed(20260719)
target_ids <- unlist(lapply(c("control", "asymptomatic", "symptomatic"), function(group)
  sample(rownames(samples)[samples$dataset == "GSE198449" & samples$clinical_group == group], 20)))
selected_ids <- c(rownames(samples)[samples$dataset == "GSE17156"], target_ids)
samples <- samples[selected_ids, , drop=FALSE]

# Produce common_expression.tsv from the two normalized matrices: intersect
# gene symbols, select the training samples, and apply the shared sample IDs.
micro <- read.delim(file.path(raw_input, "GSE17156/expression_rma.tsv"), row.names=1, check.names=FALSE)
rna <- read.delim(file.path(raw_input, "GSE198449/vst_expression.tsv"), row.names=1, check.names=FALSE)
bare_ids <- sub("^GSE[0-9]+__", "", rownames(samples))
common_genes <- intersect(rownames(micro), rownames(rna))
expression <- cbind(micro[common_genes, bare_ids[1:sum(samples$dataset == "GSE17156")], drop=FALSE],
                     rna[common_genes, bare_ids[(sum(samples$dataset == "GSE17156") + 1):length(bare_ids)], drop=FALSE])
colnames(expression) <- rownames(samples)
write.table(expression, file.path(output, "common_expression.tsv"), sep="\t", quote=FALSE, col.names=NA)
samples$coconut_control <- as.integer(samples$clinical_group != "control")

# COCONUT estimates batch correction from controls and applies the same
# correction framework to the non-control samples in each dataset.
gse_objects <- map(datasets, function(dataset) {
  ids <- rownames(samples)[samples$dataset == dataset]
  list(genes=expression[, ids, drop=FALSE], pheno=samples[ids, c("coconut_control", "clinical_group"), drop=FALSE])
}) %>% set_names(datasets)
coconut <- COCONUT(gse_objects, control.0.col="coconut_control", disease.col="clinical_group", byPlatform=FALSE, par.prior=TRUE, parallel=FALSE)

# Restore the corrected control and non-control matrices to the original sample order.
integrated <- map(datasets, function(dataset) {
  ids <- rownames(samples)[samples$dataset == dataset]
  controls <- ids[samples[ids, "coconut_control"] == 0]
  disease <- ids[samples[ids, "coconut_control"] == 1]
  corrected_control <- coconut$controlList$GSEs[[dataset]]$genes[, seq_along(controls), drop=FALSE]
  corrected_disease <- coconut$COCONUTList[[dataset]]$genes[, seq_along(disease), drop=FALSE]
  result <- cbind(corrected_control, corrected_disease); colnames(result) <- c(controls, disease); result
}) %>% reduce(cbind)
write.table(integrated, file.path(output, "coconut_expression.tsv"), sep="\t", quote=FALSE, col.names=NA)

# WGCNA retains the original 12,531-gene input and network parameters.
dat_expr <- t(integrated)
candidate_powers <- c(1:10, seq(12, 20, 2))
soft_threshold <- pickSoftThreshold(dat_expr, powerVector=candidate_powers, networkType="signed", verbose=0)
write.table(soft_threshold$fitIndices, file.path(output, "wgcna_soft_threshold.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
network <- blockwiseModules(dat_expr, power=16, networkType="signed", TOMType="signed", maxBlockSize=2000,
                            minModuleSize=500, deepSplit=2, mergeCutHeight=0.25, numericLabels=TRUE,
                            pamRespectsDendro=FALSE, saveTOMs=FALSE, verbose=0)
module_table <- data.frame(gene_symbol=colnames(dat_expr), module=labels2colors(network$colors))
write.table(module_table, file.path(output, "wgcna_modules.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Encode the nine one-hot traits used to relate module eigengenes to conditions.
condition <- ifelse(samples$dataset == "GSE198449", "S", ifelse(samples$cohort == "HRV", "H", ifelse(samples$cohort == "IFV", "I", "R")))
trait_names <- c("Control", "H_asymptomatic", "H_symptomatic", "I_asymptomatic", "I_symptomatic", "R_asymptomatic", "R_symptomatic", "S_asymptomatic", "S_symptomatic")
traits <- map_dfc(trait_names, function(x) {
  values <- if (x == "Control") samples$clinical_group == "control" else {
    p <- strsplit(x, "_")[[1]]
    condition == p[1] & samples$clinical_group == p[2]
  }
  tibble(!!x := values)
})
MEs <- orderMEs(network$MEs); correlations <- cor(MEs, traits); pvalues <- corPvalueStudent(correlations, nrow(dat_expr))
write.table(data.frame(module=rownames(correlations), correlations, check.names=FALSE), file.path(output, "wgcna_module_trait_correlations.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Intersect the WGCNA brown module with the two discovery DEG lists to define stable candidates.
target_module = "brown"
deg161918 <- read.delim(file.path(deg_input, "GSE161918/symptomatic_vs_control/deg_significant.tsv")) %>% pull(gene_symbol)
deg171110 <- read.delim(file.path(deg_input, "GSE171110/symptomatic_vs_control/deg_significant.tsv")) %>% pull(gene_symbol)
target_genes <- module_table %>% filter(module == target_module) %>% pull(gene_symbol)
intersection <- tibble(gene_symbol=intersect(target_genes, intersect(deg161918, deg171110)))
write.table(intersection, file.path(output, "target_module_deg_intersection.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Export the exact candidate matrix and manifest consumed by the next machine-learning module.
candidate_genes <- intersection %>% mutate(source_module=target_module,
  candidate_generation="sample20;minModuleSize=500;brown∩GSE161918∩GSE171110;GSE152641_excluded")
write.table(candidate_genes, file.path(output, "candidate_genes_319.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
training_expression <- t(integrated[intersection$gene_symbol, , drop=FALSE])
write.table(training_expression, file.path(output, "training_expression_319.tsv"), sep="\t", quote=FALSE, col.names=NA)
training_samples <- samples %>% rownames_to_column("sample_id") %>% mutate(binary_label=as.integer(clinical_group %in% c("asymptomatic", "symptomatic")))
write.table(training_samples, file.path(output, "training_sample_manifest.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(training_samples, file.path(output, "integrated_samples.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
