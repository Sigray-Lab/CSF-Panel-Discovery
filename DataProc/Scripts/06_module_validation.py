#!/usr/bin/env python3
"""Step 6: Module / co-abundance validation.

On D11 (N=2,720 samples):
1. Build log2 intensity matrix, filter to >= 30% completeness
2. Impute missing values (half-minimum per protein)
3. Derive co-abundance modules (correlation-based clustering)
4. Test R1 gene enrichment per module
5. Identify module-derived candidates (hubs not in R1)

Outputs:
- DerivedData/modules/module_assignments.tsv
- DerivedData/modules/module_enrichment.tsv
- DerivedData/modules/module_derived_candidates.tsv
"""

import logging
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

SCRIPTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS_DIR))

from utils.parsers import load_config, resolve_path, parse_r1_reference, identify_sample_columns

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("06_module_validation")


# ------------------------------------------------------------------ #
# Imputation
# ------------------------------------------------------------------ #

def impute_half_min(matrix: pd.DataFrame) -> pd.DataFrame:
    """Impute missing values with half the row (protein) minimum."""
    imputed = matrix.copy()
    for idx in imputed.index:
        row = imputed.loc[idx]
        row_min = row.min(skipna=True)
        if np.isnan(row_min):
            row_min = 0.0
        imputed.loc[idx] = row.fillna(row_min / 2.0)
    return imputed


# ------------------------------------------------------------------ #
# Module detection
# ------------------------------------------------------------------ #

def correlation_clustering(intensity_matrix: pd.DataFrame,
                             min_module_size: int = 20) -> pd.DataFrame:
    """Hierarchical clustering on protein-protein correlation matrix.

    Returns DataFrame with: gene_symbol, module_id, is_hub,
    module_connectivity
    """
    logger.info(f"Computing correlation matrix ({intensity_matrix.shape[0]} proteins)...")

    # Compute Pearson correlation
    corr = intensity_matrix.T.corr(method="pearson")

    # Convert to distance matrix
    dist = 1 - corr.values
    np.fill_diagonal(dist, 0)
    # Ensure symmetry and non-negative
    dist = (dist + dist.T) / 2
    dist = np.clip(dist, 0, 2)

    # Hierarchical clustering
    logger.info("Running hierarchical clustering...")
    condensed = squareform(dist, checks=False)
    Z = linkage(condensed, method="average")

    # Dynamic tree cut: try different heights to get reasonable module sizes
    best_labels = None
    best_n_modules = 0

    for height in np.arange(0.3, 0.9, 0.05):
        labels = fcluster(Z, t=height, criterion="distance")
        # Count modules with sufficient size
        module_sizes = pd.Series(labels).value_counts()
        n_good = (module_sizes >= min_module_size).sum()
        if n_good > best_n_modules and n_good <= 50:
            best_n_modules = n_good
            best_labels = labels
            best_height = height

    if best_labels is None:
        # Fallback: fixed number of clusters
        best_labels = fcluster(Z, t=20, criterion="maxclust")
        best_height = "maxclust=20"

    logger.info(f"  Best height: {best_height}, {best_n_modules} modules >= {min_module_size} proteins")

    # Build results
    genes = intensity_matrix.index.tolist()
    module_df = pd.DataFrame({
        "gene_symbol": genes,
        "module_id": best_labels,
    })

    # Compute module connectivity (mean correlation with module members)
    module_df["module_connectivity"] = 0.0
    for mod_id in module_df["module_id"].unique():
        mod_mask = module_df["module_id"] == mod_id
        mod_genes = module_df[mod_mask]["gene_symbol"].tolist()
        if len(mod_genes) < 2:
            continue

        mod_indices = [genes.index(g) for g in mod_genes if g in genes]
        for g in mod_genes:
            if g not in genes:
                continue
            g_idx = genes.index(g)
            conn = np.mean([corr.iloc[g_idx, j] for j in mod_indices if j != g_idx])
            module_df.loc[module_df["gene_symbol"] == g, "module_connectivity"] = conn

    # Hub proteins: top 10% connectivity within each module
    module_df["is_hub"] = False
    for mod_id in module_df["module_id"].unique():
        mod_mask = module_df["module_id"] == mod_id
        if mod_mask.sum() < min_module_size:
            continue
        threshold = module_df.loc[mod_mask, "module_connectivity"].quantile(0.90)
        hub_mask = mod_mask & (module_df["module_connectivity"] >= threshold)
        module_df.loc[hub_mask, "is_hub"] = True

    # Module sizes
    module_sizes = module_df.groupby("module_id").size()
    logger.info(f"  Module sizes: {dict(module_sizes[module_sizes >= min_module_size])}")

    return module_df


# ------------------------------------------------------------------ #
# Enrichment testing
# ------------------------------------------------------------------ #

def test_module_enrichment(modules: pd.DataFrame, r1_genes: set,
                             total_genes: int) -> pd.DataFrame:
    """Per-module Fisher's exact test for R1 enrichment."""
    r1_upper = {g.upper() for g in r1_genes}
    modules["in_r1"] = modules["gene_symbol"].str.upper().isin(r1_upper)

    results = []
    for mod_id, group in modules.groupby("module_id"):
        n_module = len(group)
        n_r1_in_module = group["in_r1"].sum()
        n_r1_total = modules["in_r1"].sum()
        n_not_r1_not_module = total_genes - n_module - n_r1_total + n_r1_in_module

        # Fisher's exact test (2x2 contingency table)
        table = [
            [n_r1_in_module, n_module - n_r1_in_module],
            [n_r1_total - n_r1_in_module,
             max(0, total_genes - n_module - n_r1_total + n_r1_in_module)],
        ]
        try:
            odds_ratio, p_value = stats.fisher_exact(table, alternative="greater")
        except Exception:
            odds_ratio, p_value = 1.0, 1.0

        expected = n_module * n_r1_total / max(total_genes, 1)
        fold_enrichment = n_r1_in_module / max(expected, 0.01)

        results.append({
            "module_id": mod_id,
            "n_genes": n_module,
            "n_r1": n_r1_in_module,
            "expected_r1": round(expected, 2),
            "fold_enrichment": round(fold_enrichment, 2),
            "pvalue": p_value,
        })

    enrichment = pd.DataFrame(results).sort_values("pvalue")

    # BH correction
    n_tests = len(enrichment)
    enrichment["rank"] = range(1, n_tests + 1)
    enrichment["pvalue_adj"] = enrichment["pvalue"] * n_tests / enrichment["rank"]
    enrichment["pvalue_adj"] = enrichment["pvalue_adj"].clip(upper=1.0)
    # Ensure monotonicity
    enrichment["pvalue_adj"] = enrichment["pvalue_adj"][::-1].cummin()[::-1]
    enrichment = enrichment.drop(columns=["rank"])

    n_sig = (enrichment["pvalue_adj"] < 0.05).sum()
    logger.info(f"  {n_sig} modules significantly enriched (FDR < 0.05)")

    return enrichment


def main():
    config_path = SCRIPTS_DIR / "config.yaml"
    config = load_config(config_path)
    raw_dir = (SCRIPTS_DIR / config["paths"]["raw_dir"]).resolve()
    derived_dir = (SCRIPTS_DIR / config["paths"]["derived_dir"]).resolve()
    log_dir = (SCRIPTS_DIR / config["paths"]["log_dir"]).resolve()
    mod_dir = derived_dir / "modules"
    mod_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    file_handler = logging.FileHandler(log_dir / f"06_module_validation_{timestamp}.log")
    file_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logging.getLogger().addHandler(file_handler)

    mod_config = config["modules"]

    logger.info("=" * 70)
    logger.info("Step 6: Module / Co-abundance Validation")
    logger.info("=" * 70)

    # ------------------------------------------------------------------
    # 1. Load D11 standardised data
    # ------------------------------------------------------------------
    d11_path = derived_dir / "standardised" / "D11_standardised.csv"
    if not d11_path.exists():
        raise FileNotFoundError(f"D11 standardised file not found: {d11_path}")

    logger.info("Loading D11...")
    d11 = pd.read_csv(d11_path, low_memory=False)
    logger.info(f"  D11: {d11.shape}")

    # ------------------------------------------------------------------
    # 2. Build intensity matrix
    # ------------------------------------------------------------------
    d11_config = config["datasets"]["D11"]
    sample_cols = identify_sample_columns(d11, d11_config)
    logger.info(f"  Sample columns: {len(sample_cols)}")

    # Use human_gene_symbol as index; drop duplicates (keep first)
    d11_dedup = d11.dropna(subset=["human_gene_symbol"]).drop_duplicates(
        subset=["human_gene_symbol"], keep="first"
    )
    intensity = d11_dedup.set_index("human_gene_symbol")[sample_cols].astype(np.float32)
    logger.info(f"  Intensity matrix: {intensity.shape}")

    # ------------------------------------------------------------------
    # 3. Filter by completeness
    # ------------------------------------------------------------------
    min_completeness = mod_config.get("min_completeness", 0.30)
    completeness = intensity.notna().mean(axis=1)
    intensity = intensity[completeness >= min_completeness]
    logger.info(f"  After completeness filter (>= {min_completeness}): {intensity.shape[0]} proteins")

    if intensity.shape[0] < 50:
        raise ValueError(f"Too few proteins ({intensity.shape[0]}) after completeness filtering")

    # ------------------------------------------------------------------
    # 4. Impute missing values
    # ------------------------------------------------------------------
    logger.info("Imputing missing values (half-minimum)...")
    intensity = impute_half_min(intensity)

    # Verify no NaN remaining
    n_nan = intensity.isna().sum().sum()
    if n_nan > 0:
        logger.warning(f"  {n_nan} NaN values remain after imputation, filling with 0")
        intensity = intensity.fillna(0)

    # ------------------------------------------------------------------
    # 5. Module detection
    # ------------------------------------------------------------------
    min_module_size = mod_config.get("min_module_size", 20)
    modules = correlation_clustering(intensity, min_module_size=min_module_size)
    logger.info(f"  Module assignments: {len(modules)} proteins across "
                f"{modules['module_id'].nunique()} modules")

    # ------------------------------------------------------------------
    # 6. R1 enrichment testing
    # ------------------------------------------------------------------
    r1_ref = config["references"]["R1"]
    r1_path = resolve_path(raw_dir, r1_ref["file"])
    r1 = parse_r1_reference(r1_path, r1_ref)

    enrichment = test_module_enrichment(modules, r1["genes"], len(intensity))

    # ------------------------------------------------------------------
    # 7. Identify module-derived candidates
    # ------------------------------------------------------------------
    r1_upper = {g.upper() for g in r1["genes"]}
    enriched_modules = enrichment[enrichment["pvalue_adj"] < 0.05]["module_id"].tolist()

    module_candidates = modules[
        (modules["module_id"].isin(enriched_modules))
        & (modules["is_hub"])
        & (~modules["gene_symbol"].str.upper().isin(r1_upper))
    ].copy()
    module_candidates = module_candidates.sort_values("module_connectivity", ascending=False)
    logger.info(f"  Module-derived candidates (hub + not in R1): {len(module_candidates)}")

    # ------------------------------------------------------------------
    # 8. Save outputs
    # ------------------------------------------------------------------
    modules.to_csv(mod_dir / "module_assignments.tsv", sep="\t", index=False)
    enrichment.to_csv(mod_dir / "module_enrichment.tsv", sep="\t", index=False)
    module_candidates.to_csv(mod_dir / "module_derived_candidates.tsv", sep="\t", index=False)

    # Generate module enrichment figure
    try:
        from utils.visualization import plot_module_enrichment
        fig_dir = (SCRIPTS_DIR / config["paths"]["output_dir"]).resolve() / "figures"
        fig_dir.mkdir(parents=True, exist_ok=True)
        plot_module_enrichment(enrichment, fig_dir / "module_enrichment.pdf")
    except Exception as e:
        logger.warning(f"Could not generate module enrichment plot: {e}")

    logger.info("\n[06_module_validation] DONE")


if __name__ == "__main__":
    main()
