"""Fit optimal linear scores and validate them across all cohorts."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

# All input paths are existing preprocessing artifacts; this step does not normalize expression again.
ROOT = Path(__file__).resolve().parents[2]
STEP3 = ROOT / "data/results/3_coconut_wgcna_intersection"
INPUT = ROOT / "data/analysis_inputs"
OUTPUT = ROOT / "data/results/5_optimal_linear_combination_and_auc"
RRA_INPUT = ROOT / "data/results/4_machine_learning_and_rra/rra_ranked_genes.tsv"
SIGNATURES = {"SIAH1_PALLD": ("SIAH1", "PALLD")}
POSITIVE_GROUPS = {"asymptomatic", "symptomatic"}
# One row describes one cohort; role, file type, and both input artifacts remain explicit.
COHORTS = pd.DataFrame([
    ("discovery_train_integrated", "discovery_train", "tsv", STEP3 / "coconut_expression.tsv", STEP3 / "integrated_samples.tsv"),
    ("GSE17156", "discovery_train_cohort", "h5ad", INPUT / "GSE17156/expression.h5ad", INPUT / "GSE17156/samples.tsv"),
    ("GSE198449", "substitute_discovery_train_cohort", "tsv", INPUT / "GSE198449/vst_expression.tsv", INPUT / "GSE198449/samples.tsv"),
    ("GSE152641", "discovery_test", "tsv", INPUT / "GSE152641/vst_expression.tsv", INPUT / "GSE152641/samples.tsv"),
    ("GSE161918", "discovery_test", "tsv", INPUT / "GSE161918/vst_expression.tsv", INPUT / "GSE161918/samples.tsv"),
    ("GSE171110", "discovery_test", "tsv", INPUT / "GSE171110/vst_expression.tsv", INPUT / "GSE171110/samples.tsv"),
    ("GSE157103", "validation", "h5ad", INPUT / "GSE157103/expression.h5ad", INPUT / "GSE157103/samples.tsv"),
    ("GSE166190", "validation", "h5ad", INPUT / "GSE166190/expression.h5ad", INPUT / "GSE166190/samples.tsv"),
    ("GSE201530", "validation", "h5ad", INPUT / "GSE201530/expression.h5ad", INPUT / "GSE201530/samples.tsv"),
    ("E-MTAB-10022", "validation", "h5ad", INPUT / "E-MTAB-10022/expression.h5ad", INPUT / "E-MTAB-10022/samples.tsv"),
    ("GSE38900_platform1", "validation", "h5ad", INPUT / "GSE38900/expression_platform1.h5ad", INPUT / "GSE38900/samples_platform1.tsv"),
    ("GSE38900_platform2", "validation", "h5ad", INPUT / "GSE38900/expression_platform2.h5ad", INPUT / "GSE38900/samples_platform2.tsv")
], columns=("name", "role", "kind", "expression_path", "samples_path"))


def load_values(expression_path, samples_path, genes, kind):
    """Load one cohort as samples-by-genes values plus its clinical manifest."""
    # Both loaders return the same shape: rows are sample IDs and columns are gene symbols.
    if kind == "tsv":
        expression = pd.read_csv(expression_path, sep="\t", index_col=0)
        samples = pd.read_csv(samples_path, sep="\t")
        # Step 3 writes an explicit sample_id column; other manifests use it as the first column.
        samples = samples.set_index("sample_id") if "sample_id" in samples else samples.set_index(samples.columns[0])
        # Duplicate feature rows are averaged before selecting genes in fixed order.
        values = expression.groupby(expression.index.astype(str)).mean().loc[list(genes)].T
    else:
        adata = ad.read_h5ad(expression_path)
        # H5AD features are mapped through gene_symbol, then duplicate feature columns are averaged.
        symbols = adata.var["gene_symbol"].fillna("").astype(str) if "gene_symbol" in adata.var else pd.Series(adata.var_names.astype(str), index=adata.var_names)
        columns = {}
        for gene in genes:
            indices = np.flatnonzero(symbols.to_numpy() == gene)
            selected = adata[:, indices].X
            if hasattr(selected, "toarray"):
                selected = selected.toarray()
            columns[gene] = np.asarray(selected, float).mean(axis=1)
        values, samples = pd.DataFrame(columns, index=adata.obs_names.astype(str)), adata.obs.copy()
        if "clinical_group" not in samples:
            samples = pd.read_csv(samples_path, sep="\t", index_col=0)
    samples.index = samples.index.astype(str)
    return values.loc[samples.index], samples


def labels(samples):
    """Encode the target condition as one and control/comparison groups as zero."""
    return samples["clinical_group"].astype(str).isin(POSITIVE_GROUPS).astype(int).to_numpy()


def fit_weights(values, target):
    """Estimate covariance-adjusted weights, normalize them, and orient scores by training AUC."""
    positive, negative = values.to_numpy()[target == 1], values.to_numpy()[target == 0]
    mean_difference = positive.mean(axis=0) - negative.mean(axis=0)
    covariance = np.cov(positive, rowvar=False) + np.cov(negative, rowvar=False)
    raw = np.linalg.pinv(covariance) @ mean_difference
    normalized = raw / raw.sum()  # The requested optimal linear combination is normalized to sum one.
    train_score = values.to_numpy() @ normalized
    orientation = 1 if roc_auc_score(target, train_score) >= .5 else -1
    return raw, normalized, normalized * orientation


def evaluate(name, role, values, samples, genes, weights):
    """Apply locked weights and save one score and AUC record for a cohort."""
    target = labels(samples)
    score = values.loc[samples.index, list(genes)].to_numpy() @ weights
    table = samples[["clinical_group"]].copy().assign(sample_id=samples.index, dataset=name, role=role, binary_label=target)[["sample_id", "clinical_group", "dataset", "role", "binary_label"]]
    table["linear_score"] = score
    for gene, weight in zip(genes, weights):
        table[f"{gene}_contribution"] = values.loc[samples.index, gene].to_numpy() * weight
    auc = roc_auc_score(target, score) if len(np.unique(target)) == 2 else np.nan
    result = {"dataset": name, "role": role, "n_samples": len(target), "positive_n": int(target.sum()), "negative_n": int((target == 0).sum()), "auc": auc, "status": "available" if np.isfinite(auc) else "not_estimable_one_class"}
    return result, table



# Require the previous module's RRA artifact before evaluating fixed signatures.
OUTPUT.mkdir(parents=True, exist_ok=True)
rra_genes = set(pd.read_csv(RRA_INPUT, sep="\t")["gene"])
assert all(set(genes) <= rra_genes for genes in SIGNATURES.values())
loaded = {}
# Load every cohort once so each signature uses identical samples and annotations.
for cohort in COHORTS.itertuples(index=False):
    loaded[cohort.name] = (cohort.role, *load_values(cohort.expression_path, cohort.samples_path, set().union(*SIGNATURES.values()), cohort.kind))
for signature, genes in SIGNATURES.items():
    _, train_values, train_samples = loaded["discovery_train_integrated"]
    train_values = train_values.loc[:, list(genes)]
    target = labels(train_samples)
    # Fit weights only on the integrated training cohort; all other cohorts are evaluation-only.
    raw, normalized, weights = fit_weights(train_values, target)
    folder = OUTPUT / signature
    folder.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"gene": genes, "raw_weight": raw, "normalized_weight_sum1": normalized, "evaluation_weight": weights}).to_csv(folder / "model_parameters.tsv", sep="\t", index=False)
    auc_rows, score_tables = [], []
    for name, (role, values, samples) in loaded.items():
        result, table = evaluate(name, role, values, samples, genes, weights)
        auc_rows.append(result); score_tables.append(table)
    pd.DataFrame(auc_rows).to_csv(folder / "cohort_auc.tsv", sep="\t", index=False)
    pd.concat(score_tables, ignore_index=True).to_csv(folder / "sample_linear_scores.tsv", sep="\t", index=False)
