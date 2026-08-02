# 1. DEG analysis and GO/KEGG enrichment
# Inputs are existing normalized/exported artifacts; no downloading or preprocessing is performed here.
suppressPackageStartupMessages({
  library(limma); library(DESeq2); library(clusterProfiler)
  library(org.Hs.eg.db); library(AnnotationDbi)
  library(tidyverse)
})

root <- normalizePath(".")
# All external expression and sample artifacts are staged under one teaching input tree.
input <- file.path(root, "data/analysis_inputs")
output <- file.path(root, "data/results/1_deg_and_enrichment")

# Convert symbols once so every enrichment uses the same annotation universe.
map_symbols <- function(symbols) {
  map <- AnnotationDbi::select(org.Hs.eg.db, unique(symbols), "ENTREZID", "SYMBOL")
  map %>% filter(!is.na(ENTREZID)) %>% distinct(ENTREZID) %>% pull(ENTREZID)
}

# Conduct Go and KEGG analysis
enrich_set <- function(genes, universe, folder, label) {
  ids <- map_symbols(genes); bg <- map_symbols(universe)
  walk(c("BP", "CC", "MF"), function(ontology) {
    result <- enrichGO(ids, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont=ontology, # test over-representation for one GO ontology
                       universe=bg, pAdjustMethod="BH", readable=TRUE)
    write.table(as.data.frame(result), file.path(folder, paste0("go_", tolower(ontology), "_", label, ".tsv")),
                sep="\t", row.names=FALSE, quote=FALSE)
  })
  result <- tryCatch(enrichKEGG(ids, universe=bg, organism="hsa", pAdjustMethod="BH"), # test KEGG pathway over-representation
                     error=function(e) NULL)
  write.table(if (is.null(result)) data.frame() else as.data.frame(result),
              file.path(folder, paste0("kegg_", label, ".tsv")), sep="\t", row.names=FALSE, quote=FALSE)
}

# limma handles the existing log2 RMA matrix from the microarray cohort.
run_limma <- function(expression, samples, comparison, folder) {
  keep <- samples %>% pull(clinical_group) %in% c("control", comparison)
  samples <- samples[keep, , drop=FALSE] %>% droplevels(); expression <- expression[, keep]
  samples$clinical_group <- factor(samples$clinical_group, levels=c("control", comparison))
  design <- model.matrix(~0 + clinical_group, samples) # one design column per group
  colnames(design) <- sub("clinical_group", "", colnames(design)) # expose clean contrast names
  fit <- lmFit(expression, design) # fit one linear model for every gene
  fit <- contrasts.fit(fit, makeContrasts(contrasts=paste0(comparison, "-control"), levels=design)) # define the pairwise contrast
  fit <- eBayes(fit) # moderate gene-wise variances with an empirical-Bayes prior
  result <- topTable(fit, number=Inf, sort.by="none", adjust.method="BH") # calculate logFC, p values, and BH padj
  result <- result %>% rownames_to_column("gene_symbol") %>%
    rename(log2FoldChange=logFC, pvalue=P.Value, padj=adj.P.Val) %>%
    mutate(direction=if_else(log2FoldChange > 0, "up", "down"))
  save_outputs(result, samples, comparison, folder)
}

# DESeq2 handles the existing count matrix from the RNA-seq cohort.
run_deseq2 <- function(expression, samples, comparison, folder) {
  keep <- samples %>% pull(clinical_group) %in% c("control", comparison)
  samples <- samples[keep, , drop=FALSE] %>% droplevels(); expression <- round(expression[, keep])
  samples$clinical_group <- factor(samples$clinical_group, levels=c("control", comparison))
  expression <- expression[rowSums(expression >= 10) >= 3, ]
  dds <- DESeqDataSetFromMatrix(expression, samples, ~clinical_group) # construct the DESeq2 experiment object
  dds <- DESeq(dds, quiet=TRUE) # estimate size factors, dispersion, and the count model
  result <- as.data.frame(results(dds, contrast=c("clinical_group", comparison, "control"))) # extract the requested contrast
  result <- result %>% rownames_to_column("gene_symbol") %>%
    mutate(direction=if_else(log2FoldChange > 0, "up", "down"))
  save_outputs(result, samples, comparison, folder)
}

# Save the complete DEG table, directional subsets, and three enrichment sets.
save_outputs <- function(result, samples, comparison, folder) {
  dir.create(folder, recursive=TRUE, showWarnings=FALSE)
  significant <- result %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 0)
  up <- significant %>% filter(log2FoldChange > 0); down <- significant %>% filter(log2FoldChange < 0)
  tables <- list(all=result, significant=significant, up=up, down=down)
  iwalk(tables, ~write.table(.x, file.path(folder, paste0("deg_", .y, ".tsv")), sep="\t", row.names=FALSE, quote=FALSE))
  enrich_set(significant$gene_symbol, result$gene_symbol, folder, "all")
  enrich_set(up$gene_symbol, result$gene_symbol, folder, "up")
  enrich_set(down$gene_symbol, result$gene_symbol, folder, "down")
}

# Use asymptomatic-vs-control and symptomatic-vs-control.
micro_expr <- read.delim(file.path(input, "GSE17156/expression_rma.tsv"), row.names=1, check.names=FALSE)
micro_samples <- read.delim(file.path(input, "GSE17156/samples.tsv"), row.names=1)
for (virus in unique(micro_samples$virus)) for (comparison in c("asymptomatic", "symptomatic")) {
  micro_subset <- micro_samples[micro_samples$virus == virus, , drop=FALSE]
  run_limma(as.matrix(micro_expr[, rownames(micro_subset)]), micro_subset, comparison,
            file.path(output, "GSE17156", virus, paste0(comparison, "_vs_control")))
}

# DESeq2 is applied to all three RNA-seq discovery cohorts with the same comparisons.
for (dataset in c("GSE198449", "GSE161918", "GSE171110")) {
  count_expr <- read.delim(file.path(input, dataset, "counts.tsv"), row.names=1, check.names=FALSE)
  count_samples <- read.delim(file.path(input, dataset, "samples.tsv"), row.names=1)
  comparisons <- intersect(c("asymptomatic", "symptomatic"), unique(count_samples$clinical_group))
  for (comparison in comparisons)
    run_deseq2(as.matrix(count_expr[, rownames(count_samples)]), count_samples, comparison,
               file.path(output, dataset, paste0(comparison, "_vs_control")))
}
