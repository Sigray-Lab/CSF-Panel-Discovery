#!/usr/bin/env python3
"""
Simulation: What if R1 was extended with GO/Reactome autophagy genes?

This standalone script re-scores and re-selects the core panel as if
proteins validated by Gene Ontology/Reactome as autophagy-related (but
NOT in the curated R1 list) received partial autophagy credit.

Key question: Would SEC22B (and other GO-validated autophagy genes)
enter the core panel if their curation gap were closed?

IMPORTANT: This does NOT modify any existing pipeline outputs.
All results go to  QC/sim_extended_r1/
"""

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

# ------------------------------------------------------------------ #
# Paths
# ------------------------------------------------------------------ #
SCRIPTS_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPTS_DIR.parent
QC_DIR = BASE_DIR / "QC"
DERIVED_DIR = BASE_DIR / "DerivedData"
OUTPUT_DIR = QC_DIR / "sim_extended_r1"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(OUTPUT_DIR / "sim_extended_r1.log", mode="w"),
    ],
)
logger = logging.getLogger(__name__)


# ------------------------------------------------------------------ #
# Helpers  (self-contained — no imports from pipeline utils)
# ------------------------------------------------------------------ #

def load_config():
    with open(SCRIPTS_DIR / "config.yaml") as f:
        return yaml.safe_load(f)


def load_evidence_scores() -> pd.DataFrame:
    """Load the pre-computed evidence scores from pipeline Step 3."""
    path = DERIVED_DIR / "evidence_scores" / "evidence_scores.tsv"
    return pd.read_csv(path, sep="\t")


def load_r1_validation() -> pd.DataFrame:
    """Load QC/R1_validation_vs_GO_Reactome.tsv"""
    path = QC_DIR / "R1_validation_vs_GO_Reactome.tsv"
    return pd.read_csv(path, sep="\t")


def load_original_shortlist() -> set:
    """Load the original 80-protein core panel shortlist."""
    path = BASE_DIR / "Outputs" / "core_panel_shortlist.tsv"
    df = pd.read_csv(path, sep="\t")
    return set(df["human_gene_symbol"].dropna())


def build_protein_evidence(config):
    """Rebuild protein evidence dict from standardised files (for null sim)."""
    protein_evidence = {}
    for ds_id in config["datasets"]:
        if ds_id == "EV":
            continue
        ds_path = DERIVED_DIR / "standardised" / f"{ds_id}_standardised.csv"
        if not ds_path.exists():
            continue
        df = pd.read_csv(ds_path, low_memory=False)
        for _, row in df.iterrows():
            gene = row.get("human_gene_symbol")
            if pd.isna(gene) or str(gene).strip() in ("", "nan", "None"):
                continue
            gene = str(gene).strip()
            if gene not in protein_evidence:
                protein_evidence[gene] = {
                    "is_ambiguous_group": False,
                    "orthology_ambiguous": False,
                    "likely_plasma_derived": False,
                }
            detected = (row.get("n_detected", 0) > 0 or
                        row.get("detectability_tier", "absent") != "absent")
            protein_evidence[gene][ds_id] = {
                "detected": detected,
                "n_detected": int(row.get("n_detected", 0)),
                "pct_detected": float(row.get("pct_detected", 0.0)),
                "detectability_tier": row.get("detectability_tier", "absent"),
            }
            if row.get("is_ambiguous_group", False):
                protein_evidence[gene]["is_ambiguous_group"] = True
            if row.get("orthology_ambiguous", False):
                protein_evidence[gene]["orthology_ambiguous"] = True
            if row.get("likely_plasma_derived", False):
                protein_evidence[gene]["likely_plasma_derived"] = True
    return protein_evidence


# ------------------------------------------------------------------ #
# Main simulation
# ------------------------------------------------------------------ #

def main():
    logger.info("=" * 70)
    logger.info("SIMULATION: Extended R1 (+ GO/Reactome autophagy genes)")
    logger.info("=" * 70)

    config = load_config()
    scoring_cfg = config["scoring"]
    weights = scoring_cfg["weights"]
    GO_MODIFIER = 0.70  # conservative: same as upstream_regulators tier

    # ------------------------------------------------------------------
    # 1. Load data
    # ------------------------------------------------------------------
    evidence = load_evidence_scores()
    r1_val = load_r1_validation()
    original_shortlist = load_original_shortlist()

    logger.info(f"Evidence scores: {len(evidence)} proteins")
    logger.info(f"R1 validation: {len(r1_val)} genes total")

    # Identify go_reactome_only genes
    go_only = r1_val[
        (r1_val["in_r1"] == False) &
        ((r1_val["in_go_autophagy"] == True) | (r1_val["in_reactome_autophagy"] == True))
    ]["gene_symbol"].dropna().astype(str).str.strip()
    go_only_set = set(go_only)

    # Original R1 genes
    r1_genes = set(
        r1_val[r1_val["in_r1"] == True]["gene_symbol"].dropna().astype(str).str.strip()
    )

    logger.info(f"Curated R1 genes: {len(r1_genes)}")
    logger.info(f"GO/Reactome-only genes (not in R1): {len(go_only_set)}")

    # Extended set
    extended_autophagy = r1_genes | go_only_set
    logger.info(f"Extended autophagy set: {len(extended_autophagy)}")

    # ------------------------------------------------------------------
    # 2. Re-score: give GO-only genes partial autophagy credit
    # ------------------------------------------------------------------
    evidence = evidence.copy()

    # Store original values for comparison
    evidence["original_score_autophagy"] = evidence["score_autophagy"]
    evidence["original_composite"] = evidence["composite_score"]
    evidence["original_rank"] = evidence["rank"]

    # Identify which proteins get a boost
    evidence["is_go_only"] = evidence["human_gene_symbol"].isin(go_only_set)
    evidence["in_extended_autophagy"] = evidence["human_gene_symbol"].isin(extended_autophagy)

    # Recalculate score_autophagy for go_only genes
    evidence.loc[evidence["is_go_only"], "score_autophagy"] = GO_MODIFIER

    # Recalculate composite_score for ALL proteins (formula unchanged)
    def recalc_composite(row):
        components = (
            weights["mouse_csf"] * row["score_mouse_csf"]
            + weights["human_csf"] * row["score_human_csf"]
            + weights["ev"] * row["score_ev"]
            + weights["brain_plausibility"] * row["score_brain"]
            + weights["autophagy_membership"] * row["score_autophagy"]
        )
        score = components + row["penalties"]
        return max(0.0, min(1.0, score))

    evidence["composite_score"] = evidence.apply(recalc_composite, axis=1)

    # Re-rank
    evidence = evidence.sort_values("composite_score", ascending=False).reset_index(drop=True)
    evidence["rank"] = range(len(evidence))

    # ------------------------------------------------------------------
    # 3. Apply core panel criteria (extended R1)
    # ------------------------------------------------------------------
    required_tiers = config["core_panel"].get("require_human_csf_tier", ["A", "B"])
    required_mouse_tiers = config["core_panel"].get("require_mouse_csf_tier", ["A", "B"])

    mask = pd.Series(True, index=evidence.index)

    # D11 tier A or B
    mask &= evidence["d11_tier"].isin(required_tiers)

    # Mouse CSF tier A or B
    mask &= evidence["mouse_csf_best_tier"].isin(required_mouse_tiers)

    # Extended R1 membership (original R1 OR GO/Reactome autophagy)
    mask &= evidence["in_extended_autophagy"] == True

    # Not plasma-derived
    mask &= evidence["likely_plasma_derived"] == False

    extended_panel = evidence[mask].copy()
    logger.info(f"\nExtended core panel: {len(extended_panel)} candidates "
                f"(vs original 122 with strict R1)")

    # ------------------------------------------------------------------
    # 4. Compare with original shortlist
    # ------------------------------------------------------------------
    extended_genes = set(extended_panel["human_gene_symbol"])

    new_entrants = extended_genes - original_shortlist
    dropped = original_shortlist - extended_genes
    # (dropped should be empty — extending R1 only adds, never removes)
    overlap = extended_genes & original_shortlist

    logger.info(f"Original shortlist: {len(original_shortlist)} proteins")
    logger.info(f"Extended panel total: {len(extended_panel)}")
    logger.info(f"Overlap with original: {len(overlap)}")
    logger.info(f"New entrants: {len(new_entrants)}")
    logger.info(f"Dropped from original: {len(dropped)}")

    # ------------------------------------------------------------------
    # 5. Output files
    # ------------------------------------------------------------------

    # 5a. Full extended shortlist
    out_cols = [
        "rank", "human_gene_symbol", "composite_score",
        "original_composite", "original_rank",
        "score_mouse_csf", "score_human_csf", "score_ev", "score_brain",
        "score_autophagy", "original_score_autophagy",
        "penalties", "n_mouse_csf_datasets", "mouse_csf_best_tier",
        "d11_tier", "d12_validated", "d10_detected", "d6_brain_detected",
        "ev_present", "in_r1", "is_go_only", "in_extended_autophagy",
        "likely_plasma_derived",
    ]
    # Keep only columns that exist
    out_cols = [c for c in out_cols if c in extended_panel.columns]
    extended_panel[out_cols].to_csv(
        OUTPUT_DIR / "extended_shortlist.tsv", sep="\t", index=False
    )
    logger.info(f"  Saved: extended_shortlist.tsv ({len(extended_panel)} proteins)")

    # 5b. New entrants (proteins entering due to GO extension)
    new_df = extended_panel[extended_panel["human_gene_symbol"].isin(new_entrants)]
    new_df[out_cols].to_csv(
        OUTPUT_DIR / "new_entrants.tsv", sep="\t", index=False
    )
    logger.info(f"  Saved: new_entrants.tsv ({len(new_df)} proteins)")

    # 5c. Dropped (should be empty unless rescoring changes ranks below max_size)
    if dropped:
        dropped_df = evidence[evidence["human_gene_symbol"].isin(dropped)]
        dropped_df[out_cols].to_csv(
            OUTPUT_DIR / "dropped.tsv", sep="\t", index=False
        )
        logger.info(f"  Saved: dropped.tsv ({len(dropped_df)} proteins)")
    else:
        logger.info("  No proteins dropped (expected: extension only adds)")

    # ------------------------------------------------------------------
    # 6. Null simulation — extended autophagy list
    # ------------------------------------------------------------------
    logger.info("\nRunning null simulation for extended R1...")

    protein_evidence = build_protein_evidence(config)
    mouse_ds = scoring_cfg["mouse_csf_datasets"]
    all_genes = list(protein_evidence.keys())

    sim_cfg = config["sensitivity"]["null_simulation"]
    n_iterations = sim_cfg["n_iterations"]
    seed = sim_cfg["seed"]
    rng = np.random.default_rng(seed)

    # --- Original R1 convergence ---
    original_observed = 0
    r1_upper = {g.upper() for g in r1_genes}
    for gene in all_genes:
        if gene.upper() not in r1_upper:
            continue
        ev = protein_evidence[gene]
        has_mouse = any(ev.get(ds, {}).get("detected", False) for ds in mouse_ds)
        has_human = ev.get("D11", {}).get("detected", False)
        if has_mouse and has_human:
            original_observed += 1

    # --- Extended R1 convergence ---
    extended_observed = 0
    ext_upper = {g.upper() for g in extended_autophagy}
    for gene in all_genes:
        if gene.upper() not in ext_upper:
            continue
        ev = protein_evidence[gene]
        has_mouse = any(ev.get(ds, {}).get("detected", False) for ds in mouse_ds)
        has_human = ev.get("D11", {}).get("detected", False)
        if has_mouse and has_human:
            extended_observed += 1

    logger.info(f"  Original R1 convergence: {original_observed} / {len(r1_genes)}")
    logger.info(f"  Extended convergence: {extended_observed} / {len(extended_autophagy)}")

    # --- Null: random sets of same size as extended ---
    extended_gene_set_size = len(extended_autophagy)
    null_counts_original = []
    null_counts_extended = []

    for i in range(n_iterations):
        # Sample size = original R1 (604)
        rand_orig = set(rng.choice(all_genes, size=min(len(r1_genes), len(all_genes)),
                                    replace=False))
        count_orig = 0
        for gene in rand_orig:
            ev = protein_evidence.get(gene, {})
            has_mouse = any(ev.get(ds, {}).get("detected", False) for ds in mouse_ds)
            has_human = ev.get("D11", {}).get("detected", False)
            if has_mouse and has_human:
                count_orig += 1
        null_counts_original.append(count_orig)

        # Sample size = extended set
        rand_ext = set(rng.choice(all_genes, size=min(extended_gene_set_size, len(all_genes)),
                                   replace=False))
        count_ext = 0
        for gene in rand_ext:
            ev = protein_evidence.get(gene, {})
            has_mouse = any(ev.get(ds, {}).get("detected", False) for ds in mouse_ds)
            has_human = ev.get("D11", {}).get("detected", False)
            if has_mouse and has_human:
                count_ext += 1
        null_counts_extended.append(count_ext)

    # Statistics
    orig_null_mean = np.mean(null_counts_original)
    orig_null_std = np.std(null_counts_original)
    orig_p = sum(1 for c in null_counts_original if c >= original_observed) / n_iterations
    orig_fold = original_observed / max(orig_null_mean, 0.01)

    ext_null_mean = np.mean(null_counts_extended)
    ext_null_std = np.std(null_counts_extended)
    ext_p = sum(1 for c in null_counts_extended if c >= extended_observed) / n_iterations
    ext_fold = extended_observed / max(ext_null_mean, 0.01)

    logger.info(f"\n  Original R1 null sim:")
    logger.info(f"    Observed: {original_observed}, Null mean: {orig_null_mean:.1f} ± {orig_null_std:.1f}")
    logger.info(f"    p = {orig_p:.4f}, fold enrichment = {orig_fold:.2f}x")
    logger.info(f"  Extended R1 null sim:")
    logger.info(f"    Observed: {extended_observed}, Null mean: {ext_null_mean:.1f} ± {ext_null_std:.1f}")
    logger.info(f"    p = {ext_p:.4f}, fold enrichment = {ext_fold:.2f}x")

    # 6b. Save null simulation comparison
    null_comparison = pd.DataFrame({
        "metric": [
            "gene_set_size",
            "observed_convergent",
            "null_mean",
            "null_std",
            "p_value",
            "fold_enrichment",
        ],
        "original_r1": [
            len(r1_genes),
            original_observed,
            round(orig_null_mean, 2),
            round(orig_null_std, 2),
            round(orig_p, 4),
            round(orig_fold, 2),
        ],
        "extended_r1_with_go": [
            len(extended_autophagy),
            extended_observed,
            round(ext_null_mean, 2),
            round(ext_null_std, 2),
            round(ext_p, 4),
            round(ext_fold, 2),
        ],
    })
    null_comparison.to_csv(
        OUTPUT_DIR / "null_simulation_comparison.tsv", sep="\t", index=False
    )
    logger.info(f"  Saved: null_simulation_comparison.tsv")

    # ------------------------------------------------------------------
    # 7. Summary text
    # ------------------------------------------------------------------
    sec22b_row = extended_panel[extended_panel["human_gene_symbol"] == "SEC22B"]
    sec22b_status = "IN EXTENDED PANEL" if len(sec22b_row) > 0 else "NOT IN EXTENDED PANEL"
    sec22b_detail = ""
    if len(sec22b_row) > 0:
        r = sec22b_row.iloc[0]
        sec22b_detail = (
            f"  Rank: {r['rank']} (was {r['original_rank']})\n"
            f"  Score: {r['composite_score']:.4f} (was {r['original_composite']:.4f})\n"
            f"  Autophagy score: {r['score_autophagy']:.2f} (was {r['original_score_autophagy']:.2f})\n"
            f"  D11 tier: {r['d11_tier']}, Mouse CSF tier: {r['mouse_csf_best_tier']}\n"
            f"  Detected in {r['n_mouse_csf_datasets']} mouse CSF datasets"
        )

    # Check top 80 of extended panel
    top80 = extended_panel.head(80)
    top80_genes = set(top80["human_gene_symbol"])
    sec22b_in_top80 = "SEC22B" in top80_genes
    new_in_top80 = top80_genes - original_shortlist
    lost_from_top80 = original_shortlist - set(extended_panel.head(80)["human_gene_symbol"])

    summary = f"""SIMULATION SUMMARY: Extended R1 with GO/Reactome Autophagy Genes
======================================================================

RATIONALE
---------
The curated R1 reference list (604 genes) defines autophagy/lysosome pathway
membership and contributes 25% of the composite score. Proteins NOT in R1 get
score_autophagy = 0, even if they are recognized autophagy genes in Gene Ontology
or Reactome. This creates a hard boundary: genuine autophagy proteins missed by
the curators (like SEC22B) cannot enter the core panel.

This simulation asks: what happens if we extend R1 with GO/Reactome-validated
autophagy genes and give them conservative partial credit (modifier = {GO_MODIFIER})?

PARAMETERS
----------
- GO/Reactome-only autophagy genes identified: {len(go_only_set)}
- Autophagy modifier for GO-only genes: {GO_MODIFIER} (same as upstream_regulators tier)
- Original R1 genes: {len(r1_genes)}
- Extended autophagy gene set: {len(extended_autophagy)}
- All other scoring weights and penalties: UNCHANGED

RESULTS
-------
Original core panel: {len(original_shortlist)} proteins (with max_size=80 cap)
Extended core panel (all qualifying): {len(extended_panel)} proteins

  - Overlap with original: {len(overlap)}
  - New entrants: {len(new_entrants)}
  - Dropped from original: {len(dropped)}

Top 80 of extended panel:
  - New in top 80: {len(new_in_top80)} proteins: {', '.join(sorted(new_in_top80)) if new_in_top80 else 'none'}
  - Lost from top 80: {len(lost_from_top80)} proteins: {', '.join(sorted(lost_from_top80)) if lost_from_top80 else 'none'}

SEC22B STATUS: {sec22b_status}
{sec22b_detail}
SEC22B in top 80: {'YES' if sec22b_in_top80 else 'NO'}

NULL SIMULATION BENCHMARK
-------------------------
                          Original R1      Extended R1+GO
Gene set size:            {len(r1_genes):<20}{len(extended_autophagy)}
Observed convergence:     {original_observed:<20}{extended_observed}
Null mean (1000 iter):    {orig_null_mean:<20.1f}{ext_null_mean:.1f}
Null std:                 {orig_null_std:<20.1f}{ext_null_std:.1f}
p-value:                  {orig_p:<20.4f}{ext_p:.4f}
Fold enrichment:          {orig_fold:<20.2f}{ext_fold:.2f}

INTERPRETATION
--------------
{'The extended gene set maintains significant enrichment above random chance.' if ext_p < 0.05 else 'WARNING: Extended gene set enrichment is not significant at p < 0.05.'}
{'The fold enrichment remains comparable to the original R1.' if abs(ext_fold - orig_fold) / orig_fold < 0.2 else 'NOTE: Fold enrichment changed substantially — the GO extension may be diluting specificity.'}

NEW ENTRANTS (proteins entering panel due to GO extension):
{chr(10).join(f'  {g}' for g in sorted(new_entrants)) if new_entrants else '  (none)'}

FILES GENERATED
---------------
  extended_shortlist.tsv        Full extended panel ({len(extended_panel)} proteins)
  new_entrants.tsv              Proteins new to panel ({len(new_entrants)})
  null_simulation_comparison.tsv  Side-by-side null simulation stats
  summary.txt                   This file

NOTE: No existing pipeline outputs were modified.
"""

    with open(OUTPUT_DIR / "summary.txt", "w") as f:
        f.write(summary)
    logger.info(f"\n  Saved: summary.txt")

    # Print summary to stdout
    print("\n" + summary)


if __name__ == "__main__":
    main()
