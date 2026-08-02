"""Five-model feature selection from the frozen 319-gene training matrix.

The script intentionally uses fixed, lightweight model settings for teaching:
cross-validation estimates stability and AUC, while RRA combines ranks only.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier
from robustrankaggregpy.aggregate_ranks import aggregate_ranks
from scipy.stats import false_discovery_control

SEED = 20260719
TOP_K = 32
ROOT = Path(__file__).resolve().parents[2]
# Step 3 exports the frozen candidate set, training matrix, and sample manifest.
INPUT = ROOT / "data/results/3_coconut_wgcna_intersection"
OUTPUT = ROOT / "data/results/4_machine_learning_and_rra"

# The rows are samples and the columns are the 319 frozen candidate genes.
expression = pd.read_csv(INPUT / "training_expression_319.tsv", sep="\t", index_col=0)
manifest = pd.read_csv(INPUT / "training_sample_manifest.tsv", sep="\t")
candidates = pd.read_csv(INPUT / "candidate_genes_319.tsv", sep="\t")["gene_symbol"].tolist()

# x is the numeric model matrix; y is the binary target: positive condition versus all controls.
# The manifest lookup enforces exactly the same sample order for expression and labels.
x = expression.loc[:, candidates].to_numpy(float)
y = manifest.set_index("sample_id").loc[expression.index, "binary_label"].to_numpy(int)
candidate_count = len(candidates)
OUTPUT.mkdir(parents=True, exist_ok=True)

# Tree models rank genes by importance; linear models rank genes by absolute coefficients.
# StandardScaler is inside the linear pipelines, so scaling is refitted independently in each fold.
models = {
    "random_forest": RandomForestClassifier(n_estimators=300, max_features="sqrt", class_weight="balanced", random_state=SEED, n_jobs=1),
    "xgboost": XGBClassifier(n_estimators=200, max_depth=3, learning_rate=.1, subsample=.8, colsample_bytree=.7, eval_metric="auc", tree_method="hist", random_state=SEED, n_jobs=1, verbosity=0),
    "gbm": GradientBoostingClassifier(n_estimators=300, max_depth=2, learning_rate=.05, min_samples_leaf=5, random_state=SEED),
    "lasso": make_pipeline(StandardScaler(), LogisticRegression(solver="liblinear", l1_ratio=1, C=1, class_weight="balanced", max_iter=5000, random_state=SEED)),
    "elastic_net": make_pipeline(StandardScaler(), LogisticRegression(solver="saga", l1_ratio=.5, C=1, class_weight="balanced", max_iter=5000, random_state=SEED)),
}
METHODS = tuple(models)
def importance(model, name, valid_x, valid_y):
    """Return one importance value per candidate gene, in the original column order."""
    if name == "random_forest":
        # Permuting a gene and measuring the AUC drop reflects its predictive contribution.
        return permutation_importance(model, valid_x, valid_y, scoring="roc_auc", n_repeats=3, random_state=SEED).importances_mean
    fitted = model[-1] if hasattr(model, "steps") else model
    if name == "xgboost":
        # XGBoost calls columns f0, f1, ... when it receives a NumPy matrix; gain is the primary rank.
        gain = fitted.get_booster().get_score(importance_type="gain")
        return np.array([gain.get(f"f{i}", 0.0) for i in range(candidate_count)])
    if name == "gbm":
        # Gradient boosting's relative influence is exposed as feature_importances_.
        return fitted.feature_importances_
    # LASSO and Elastic Net use signed coefficients; sign is retained for direction reporting.
    return fitted.coef_[0]
def selected_indices(values, name):
    """Return indices selected in one fold: top 32 for trees, non-zero coefficients for regressions."""
    order = np.argsort(-np.nan_to_num(np.abs(values), nan=-np.inf))
    if name in METHODS[:3]:
        return order[:TOP_K]
    return order[np.abs(values[order]) > 1e-12][:TOP_K]


# The three-fold stratified split keeps both classes represented and gives a short stability estimate.
# Each model is cloned before fitting, preventing one fold's fitted state from leaking into another.
folds = list(StratifiedKFold(3, shuffle=True, random_state=SEED).split(x, y))
selected_frequency = {name: np.zeros(candidate_count) for name in METHODS}
performance, selected_tables = [], []

# Fit each method on every training fold; validation AUC is computed only on that fold's held-out samples.
for name, prototype in models.items():
    fold_importances, fold_auc = [], []
    for fold_id, (train, valid) in enumerate(folds, 1):
        model = clone(prototype)
        model.fit(x[train], y[train])
        scores = model.predict_proba(x[valid])[:, 1]
        fold_auc.append(roc_auc_score(y[valid], scores))
        values = importance(model, name, x[valid], y[valid])
        fold_importances.append(values)
        selected_frequency[name][selected_indices(values, name)] += 1
        performance.append({"method": name, "fold": fold_id, "roc_auc": fold_auc[-1]})

    # Average fold-specific importance before producing the method-level complete ranking.
    values = np.nanmean(fold_importances, axis=0)
    order = np.argsort(-np.nan_to_num(np.abs(values), nan=-np.inf))
    coefficient = values[order] if name in METHODS[3:] else np.full(candidate_count, np.nan)
    ranking = pd.DataFrame({"gene": np.array(candidates)[order], "method": name,
                            "rank": np.arange(1, candidate_count + 1), "importance": np.abs(values[order]),
                            "coefficient": coefficient, "selection_frequency": selected_frequency[name][order] / len(folds),
                            "validation_auc": np.mean(fold_auc)})
    ranking["direction"] = np.select([ranking.coefficient > 0, ranking.coefficient < 0], ["positive", "negative"], default="not_applicable")
    ranking["selected"] = ranking["rank"].le(TOP_K) & ((name in METHODS[:3]) | ranking["coefficient"].abs().gt(1e-12))
    ranking.to_csv(OUTPUT / f"{name}_ranked_genes.tsv", sep="\t", index=False)
    selected_tables.append(ranking.loc[ranking["selected"]])

# RRA evaluates concordant ranks, not incomparable raw importance values from different model families.
rank_tables = {name: pd.read_csv(OUTPUT / f"{name}_ranked_genes.tsv", sep="\t").set_index("gene") for name in METHODS}
# Pass only the primary top-32 lists, while telling the package that the universe contains 319 genes.
# top_cutoff corrects the null distribution for partial lists; omitted genes receive no RRA evidence.
top_lists = [table.sort_values("rank").head(TOP_K).index.tolist() for table in rank_tables.values()]
top_fraction = min(TOP_K, candidate_count) / candidate_count
rra_scores = aggregate_ranks(rank_lists=top_lists, ranked_elements=candidate_count,
                             top_cutoff=np.repeat(top_fraction, len(METHODS)))
rra = []
for gene in candidates:
    ranks = np.array([rank_tables[name].loc[gene, "rank"] for name in METHODS])
    pvalue = float(rra_scores.get(gene, 1.0))
    rra.append({"gene": gene, **{f"{name}_rank": rank_tables[name].loc[gene, "rank"] for name in METHODS},
                "rra_score": pvalue, "raw_pvalue": pvalue, "method_count": int((ranks <= TOP_K).sum()),
                "mean_rank": ranks.mean(), "best_rank": ranks.min(),
                "selection_frequency": np.mean([rank_tables[name].loc[gene, "selection_frequency"] for name in METHODS])})
rra = pd.DataFrame(rra).sort_values(["raw_pvalue", "best_rank", "gene"]).reset_index(drop=True)

# SciPy applies Benjamini-Hochberg correction and preserves the input gene order.
rra["bh_qvalue"] = false_discovery_control(rra["raw_pvalue"].to_numpy(), method="bh")
rra.to_csv(OUTPUT / "rra_ranked_genes.tsv", sep="\t", index=False)
rra[rra.bh_qvalue < .05].to_csv(OUTPUT / "rra_selected_genes_q05.tsv", sep="\t", index=False)
for k in (10, 20, 32, 50):
    rra.head(k).to_csv(OUTPUT / f"rra_top{k}.tsv", sep="\t", index=False)
pd.concat(selected_tables).to_csv(OUTPUT / "method_selected_genes.tsv", sep="\t", index=False)
pd.DataFrame(performance).to_csv(OUTPUT / "method_cv_performance.tsv", sep="\t", index=False)
pd.DataFrame({"gene": candidates, **{name: selected_frequency[name] / len(folds) for name in METHODS}}).to_csv(OUTPUT / "gene_selection_frequency.tsv", sep="\t", index=False)
