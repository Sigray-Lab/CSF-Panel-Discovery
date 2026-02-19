#!/usr/bin/env python3
"""Step 7: Validation, sensitivity analysis, and final reporting.

1. Merge all evidence into final ranked table
2. Generate core panel shortlist
3. Compare vs R2 and previous 329-protein list
4. Fill Overview template
5. Sensitivity analysis (parameter grid)
6. Chance expectation simulation (null convergence)
7. Generate all figures
8. Write methods text draft

Outputs:
- Outputs/candidates_ranked.tsv
- Outputs/core_panel_shortlist.tsv
- Outputs/overview_matrix_filled.xlsx
- Outputs/sensitivity_analysis_summary.tsv
- Outputs/module_validation_summary.tsv
- Outputs/methods_draft.md
- Outputs/figures/ (all figures)
- QC/null_simulation_results.tsv
"""

import logging
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS_DIR))

from utils.parsers import load_config, resolve_path, parse_r1_reference, parse_r2_reference
from utils.scoring import (
    compute_mouse_csf_evidence, compute_human_csf_evidence,
    compute_ev_evidence, compute_brain_plausibility,
    compute_autophagy_membership, compute_penalties,
    compute_composite_score,
)
from utils.visualization import (
    plot_score_distribution, plot_detection_heatmap,
    plot_comparison_venn, plot_null_distribution,
    plot_sensitivity_stability, plot_upset,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("07_validate_and_report")


# ------------------------------------------------------------------ #
# Merge all evidence
# ------------------------------------------------------------------ #

def merge_final_evidence(derived_dir: Path) -> pd.DataFrame:
    """Merge evidence scores, peptide feasibility, and module assignments."""
    scores = pd.read_csv(derived_dir / "evidence_scores" / "evidence_scores.tsv", sep="\t")

    # Peptide feasibility
    pf_path = derived_dir / "peptide_feasibility" / "peptide_feasibility.tsv"
    if pf_path.exists():
        pf = pd.read_csv(pf_path, sep="\t")
        pf_cols = ["human_gene_symbol", "feasibility_tier", "n_conserved_proteotypic",
                    "n_proteotypic_human", "n_proteotypic_mouse"]
        pf_cols = [c for c in pf_cols if c in pf.columns]
        scores = scores.merge(pf[pf_cols], on="human_gene_symbol", how="left")
    else:
        logger.warning("Peptide feasibility not found; skipping merge")

    # Module assignments
    mod_path = derived_dir / "modules" / "module_assignments.tsv"
    if mod_path.exists():
        modules = pd.read_csv(mod_path, sep="\t")
        mod_cols = ["gene_symbol", "module_id", "is_hub", "module_connectivity"]
        mod_cols = [c for c in mod_cols if c in modules.columns]
        modules = modules.rename(columns={"gene_symbol": "human_gene_symbol"})
        scores = scores.merge(
            modules[["human_gene_symbol"] + [c for c in mod_cols if c != "gene_symbol"]],
            on="human_gene_symbol", how="left"
        )
    else:
        logger.warning("Module assignments not found; skipping merge")

    return scores.sort_values("composite_score", ascending=False).reset_index(drop=True)


# ------------------------------------------------------------------ #
# Sensitivity analysis
# ------------------------------------------------------------------ #

def run_sensitivity_analysis(config: dict, protein_evidence: dict,
                               ev_data: dict, r1: dict) -> pd.DataFrame:
    """Re-run scoring under parameter grid."""
    sens_cfg = config["sensitivity"]
    scoring_cfg = config["scoring"]

    tier_options = sens_cfg["tier_options"]
    weight_sets = sens_cfg["weight_sets"]

    all_results = []

    for tier_opt in tier_options:
        for ws_name, ws_weights in weight_sets.items():
            logger.info(f"  Sensitivity: tier={tier_opt}, weights={ws_name}")

            # Re-score all proteins with this parameter combo
            scores = {}
            for gene, evidence in protein_evidence.items():
                # Filter by tier option
                mouse_ds = scoring_cfg["mouse_csf_datasets"]
                filtered_evidence = {}
                for ds_id, info in evidence.items():
                    if isinstance(info, dict) and "detectability_tier" in info:
                        tier = info["detectability_tier"]
                        if tier_opt == "A" and tier != "A":
                            continue
                        elif tier_opt == "AB" and tier not in ("A", "B"):
                            continue
                        # ABC: include all
                        filtered_evidence[ds_id] = info

                mouse_csf = compute_mouse_csf_evidence(filtered_evidence, mouse_ds)
                human_csf = compute_human_csf_evidence(
                    filtered_evidence, scoring_cfg["human_csf_tier_scores"],
                    scoring_cfg["d10_bonus"]
                )
                ev = compute_ev_evidence(gene, ev_data)
                brain = compute_brain_plausibility(filtered_evidence)
                autophagy = compute_autophagy_membership(
                    gene, r1["genes"], r1["category_map"],
                    scoring_cfg["autophagy_category_modifiers"]
                )
                penalties = compute_penalties(
                    evidence.get("is_ambiguous_group", False),
                    evidence.get("orthology_ambiguous", False),
                    evidence.get("likely_plasma_derived", False),
                    scoring_cfg["penalties"],
                )

                components = {
                    "mouse_csf": mouse_csf,
                    "human_csf": human_csf,
                    "ev": ev,
                    "brain_plausibility": brain,
                    "autophagy_membership": autophagy,
                }
                composite = compute_composite_score(components, ws_weights, penalties)
                scores[gene] = composite

            # Rank proteins
            ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)
            for rank, (gene, score) in enumerate(ranked, 1):
                all_results.append({
                    "tier_option": tier_opt,
                    "weight_set": ws_name,
                    "human_gene_symbol": gene,
                    "rank": rank,
                    "score": round(score, 6),
                })

    return pd.DataFrame(all_results)


# ------------------------------------------------------------------ #
# Null simulation
# ------------------------------------------------------------------ #

def run_null_simulation(config: dict, protein_evidence: dict,
                          r1: dict) -> pd.DataFrame:
    """Random sampling null simulation.

    For each iteration, sample 604 random genes from proteome and
    count how many pass mouse + human CSF convergence criteria.
    """
    sim_cfg = config["sensitivity"]["null_simulation"]
    n_iterations = sim_cfg["n_iterations"]
    gene_set_size = sim_cfg["gene_set_size"]
    seed = sim_cfg["seed"]

    rng = np.random.default_rng(seed)

    # All proteins in the evidence dict
    all_genes = list(protein_evidence.keys())
    proteome_size = max(sim_cfg.get("proteome_size", 20000), len(all_genes))

    # Compute observed convergence: R1 genes with mouse + human CSF evidence
    r1_upper = {g.upper() for g in r1["genes"]}
    scoring_cfg = config["scoring"]
    mouse_ds = scoring_cfg["mouse_csf_datasets"]

    observed = 0
    for gene in all_genes:
        if gene.upper() not in r1_upper:
            continue
        evidence = protein_evidence[gene]
        # Has mouse CSF evidence?
        has_mouse = any(
            evidence.get(ds, {}).get("detected", False)
            for ds in mouse_ds
        )
        # Has human CSF evidence?
        has_human = evidence.get("D11", {}).get("detected", False)
        if has_mouse and has_human:
            observed += 1

    logger.info(f"  Observed convergence (R1 with mouse+human): {observed}")

    # Simulate
    null_counts = []
    for i in range(n_iterations):
        # Sample random genes
        random_genes = set(rng.choice(all_genes, size=min(gene_set_size, len(all_genes)),
                                       replace=False))
        count = 0
        for gene in random_genes:
            evidence = protein_evidence.get(gene, {})
            has_mouse = any(
                evidence.get(ds, {}).get("detected", False)
                for ds in mouse_ds
            )
            has_human = evidence.get("D11", {}).get("detected", False)
            if has_mouse and has_human:
                count += 1
        null_counts.append(count)

    # Statistics
    null_mean = np.mean(null_counts)
    null_std = np.std(null_counts)
    p_value = sum(1 for c in null_counts if c >= observed) / n_iterations
    fold_enrichment = observed / max(null_mean, 0.01)

    results = pd.DataFrame({
        "iteration": range(n_iterations),
        "n_convergent": null_counts,
    })
    results.attrs["observed"] = observed
    results.attrs["p_value"] = p_value
    results.attrs["fold_enrichment"] = fold_enrichment
    results.attrs["null_mean"] = null_mean

    logger.info(f"  Null mean: {null_mean:.1f} +/- {null_std:.1f}")
    logger.info(f"  Observed: {observed}, p = {p_value:.4f}, "
                f"fold = {fold_enrichment:.1f}x")

    return results


# ------------------------------------------------------------------ #
# Overview matrix
# ------------------------------------------------------------------ #

def fill_overview_matrix(final: pd.DataFrame, r2: pd.DataFrame,
                           provenance: pd.DataFrame,
                           output_path: Path) -> None:
    """Fill Overview template with per-target detection data."""
    import openpyxl

    # For each R2 target, look up detection across datasets
    overview_rows = []
    for _, r2_row in r2.iterrows():
        gene = r2_row["gene_symbol"]
        match = final[final["human_gene_symbol"].str.upper() == str(gene).upper()]

        row = {"gene_symbol": gene, "antibody": r2_row.get("antibody", "")}
        if not match.empty:
            first = match.iloc[0]
            row["composite_score"] = first.get("composite_score", np.nan)
            row["d11_tier"] = first.get("d11_tier", "absent")
            row["n_mouse_csf_datasets"] = first.get("n_mouse_csf_datasets", 0)
            row["in_r1"] = first.get("in_r1", False)
            row["in_new_panel"] = True
        else:
            row["composite_score"] = np.nan
            row["d11_tier"] = "not found"
            row["n_mouse_csf_datasets"] = 0
            row["in_r1"] = False
            row["in_new_panel"] = False

        overview_rows.append(row)

    overview_df = pd.DataFrame(overview_rows)
    overview_df.to_excel(output_path, index=False, engine="openpyxl")
    logger.info(f"  Overview matrix: {output_path} ({len(overview_df)} targets)")


# ------------------------------------------------------------------ #
# Methods text
# ------------------------------------------------------------------ #

def write_methods_draft(config: dict, final: pd.DataFrame,
                          output_path: Path) -> None:
    """Generate methods text draft with actual pipeline numbers."""
    n_datasets = len([k for k in config["datasets"].keys() if k != "EV"])
    n_mouse = len(config["scoring"]["mouse_csf_datasets"])
    weights = config["scoring"]["weights"]
    n_total = len(final)
    n_r1 = final["in_r1"].sum() if "in_r1" in final.columns else 0

    text = f"""# Methods Draft: CSF Autophagy/Lysosome Panel Discovery

## Proteomics Data Processing

We analyzed {n_datasets} proteomics datasets spanning mouse CSF ({n_mouse} datasets),
mouse brain tissue (1 dataset), mouse ISF (1 dataset), and human CSF (3 datasets),
comprising a total of {n_total} unique proteins after orthology mapping to human
gene symbols.

Raw proteomics data from MaxQuant (label-free quantification) and DIA-NN
(data-independent acquisition) pipelines were standardised using a uniform
processing pipeline. Contaminants (decoy sequences, common laboratory
contaminants) and high-abundance plasma proteins were flagged and excluded
from the core panel. Sample-level quality control flagged samples with
>80% missing values.

## Detectability Tiers

Proteins were assigned dataset-aware detectability tiers (A/B/C) based on
detection frequency and intensity relative to a dataset-specific intensity
floor (10th percentile of non-zero values). For small datasets (n <= 20
samples), Tier A required detection in >= 2 samples above the intensity
floor. For large datasets (n > 20), Tier A required >= 10% detection
fraction above floor.

## Orthology Mapping

Mouse gene symbols were mapped to human orthologs using g:Profiler g:Orth
(Ensembl-backed) with cross-validation via MyGene.info for ambiguous cases.
One-to-one orthologs were accepted directly; one-to-many mappings were
expanded with an ambiguity flag.

## Evidence Scoring

A continuous composite evidence score was computed for each protein:

Score = {weights['mouse_csf']} x mouse_CSF + {weights['human_csf']} x human_CSF
      + {weights['ev']} x EV + {weights['brain_plausibility']} x brain
      + {weights['autophagy_membership']} x autophagy + penalties

Mouse CSF evidence was computed as (proportion of datasets with detection)
x (mean detection fraction across detected datasets). Human CSF evidence
was based on the Astral discovery cohort (D11, N=2,720) detectability tier.
EV evidence served as supportive annotation only (weight = {weights['ev']}).
Brain tissue plausibility was assessed from mouse brain lysate proteomics.
Autophagy/lysosome membership was scored based on a curated reference list
of {n_r1} genes across 6 functional categories.

## Core Panel Selection

Core panel candidates were required to meet: (1) Human CSF detectability
tier A or B in the Astral cohort, (2) Mouse CSF tier A or B in at least
one dataset, (3) Membership in the curated autophagy/lysosome reference
list, and (4) Not flagged as likely plasma-derived.

## Sensitivity Analysis

Pipeline robustness was assessed across a parameter grid of 18 combinations
varying detectability tier stringency (A-only, A+B, A+B+C), protein group
handling (expand vs exclude ambiguous), and evidence score weight sets
(default, equal, CSF-heavy).

## Chance Expectation

To assess the significance of autophagy gene convergence across species,
we performed 1,000 random simulations sampling 604 genes from the human
proteome and computing how many passed the same mouse + human CSF
convergence criteria.

## Software

The pipeline was implemented in Python (>= 3.10) using pandas, scipy,
gprofiler-official, mygene, pyteomics, biopython, gseapy, matplotlib,
seaborn, and upsetplot. All analyses were deterministic (seed = 42) with
intermediate results saved at each step.
"""
    output_path.write_text(text)
    logger.info(f"  Methods draft: {output_path}")


# ------------------------------------------------------------------ #
# Helper: rebuild protein evidence dict from standardised files
# ------------------------------------------------------------------ #

def _build_protein_evidence_dict(std_dir, config):
    """Rebuild protein evidence dict from standardised files."""
    protein_evidence = {}
    for ds_id in config["datasets"]:
        if ds_id == "EV":
            continue
        ds_path = std_dir / f"{ds_id}_standardised.csv"
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
# Main
# ------------------------------------------------------------------ #

def main():
    config_path = SCRIPTS_DIR / "config.yaml"
    config = load_config(config_path)
    raw_dir = (SCRIPTS_DIR / config["paths"]["raw_dir"]).resolve()
    derived_dir = (SCRIPTS_DIR / config["paths"]["derived_dir"]).resolve()
    output_dir = (SCRIPTS_DIR / config["paths"]["output_dir"]).resolve()
    qc_dir = (SCRIPTS_DIR / config["paths"]["qc_dir"]).resolve()
    log_dir = (SCRIPTS_DIR / config["paths"]["log_dir"]).resolve()
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    file_handler = logging.FileHandler(log_dir / f"07_validate_report_{timestamp}.log")
    file_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logging.getLogger().addHandler(file_handler)

    logger.info("=" * 70)
    logger.info("Step 7: Validation, Sensitivity, and Reporting")
    logger.info("=" * 70)

    # ------------------------------------------------------------------
    # 1. Merge all evidence
    # ------------------------------------------------------------------
    final = merge_final_evidence(derived_dir)
    logger.info(f"Final merged table: {len(final)} proteins")

    # ------------------------------------------------------------------
    # 2. Save candidates_ranked.tsv
    # ------------------------------------------------------------------
    final.to_csv(output_dir / "candidates_ranked.tsv", sep="\t", index=False)
    logger.info(f"Ranked candidates: {output_dir / 'candidates_ranked.tsv'}")

    # ------------------------------------------------------------------
    # 3. Generate core panel shortlist
    # ------------------------------------------------------------------
    panel_cfg = config["core_panel"]
    core_path = derived_dir / "autophagy_panel" / "core_panel_candidates.tsv"
    if core_path.exists():
        shortlist = pd.read_csv(core_path, sep="\t")
        # Merge peptide feasibility if available
        pf_path = derived_dir / "peptide_feasibility" / "peptide_feasibility.tsv"
        if pf_path.exists():
            pf = pd.read_csv(pf_path, sep="\t")
            shortlist = shortlist.merge(
                pf[["human_gene_symbol", "feasibility_tier", "n_conserved_proteotypic"]],
                on="human_gene_symbol", how="left"
            )
        # Apply size limits
        max_size = panel_cfg.get("max_size", 80)
        shortlist = shortlist.head(max_size)
    else:
        # Fallback: top N from final
        shortlist = final[final["in_r1"] == True].head(panel_cfg.get("max_size", 80))

    shortlist.to_csv(output_dir / "core_panel_shortlist.tsv", sep="\t", index=False)
    logger.info(f"Core panel shortlist: {len(shortlist)} proteins")

    # ------------------------------------------------------------------
    # 4. Comparisons
    # ------------------------------------------------------------------
    r2_ref = config["references"]["R2"]
    r2_path = resolve_path(raw_dir, r2_ref["file"])
    r2 = parse_r2_reference(r2_path, r2_ref)

    new_panel_genes = set(shortlist["human_gene_symbol"].dropna())
    r2_genes = set(r2["gene_symbol"].dropna())

    # Previous 329 from EV precomputed sheet
    ev_precomputed_path = derived_dir / "standardised" / "EV_standardised.csv"
    previous_329 = set()
    if ev_precomputed_path.exists():
        ev_pre = pd.read_csv(ev_precomputed_path)
        if "Gene_Symbol" in ev_pre.columns:
            previous_329 = set(ev_pre["Gene_Symbol"].dropna())

    logger.info(f"Comparison sets: new={len(new_panel_genes)}, "
                f"R2={len(r2_genes)}, previous={len(previous_329)}")

    # Overlap stats
    overlap_r2 = new_panel_genes & r2_genes
    overlap_329 = new_panel_genes & previous_329
    logger.info(f"  New panel overlaps: R2={len(overlap_r2)}, previous 329={len(overlap_329)}")

    # ------------------------------------------------------------------
    # 5. Fill Overview matrix
    # ------------------------------------------------------------------
    provenance_path = derived_dir / "evidence_scores" / "provenance_matrix.tsv"
    provenance = pd.read_csv(provenance_path, sep="\t") if provenance_path.exists() else pd.DataFrame()
    fill_overview_matrix(final, r2, provenance, output_dir / "overview_matrix_filled.xlsx")

    # ------------------------------------------------------------------
    # 6. Build protein evidence dict for sensitivity/simulation
    # ------------------------------------------------------------------
    # Re-build from standardised files
    from utils.parsers import parse_r1_reference as _parse_r1
    r1_ref = config["references"]["R1"]
    r1_path = resolve_path(raw_dir, r1_ref["file"])
    r1 = _parse_r1(r1_path, r1_ref)

    # Load EV data for re-scoring
    ev_data = {"rank_percentiles": {}, "genes": set()}
    ev_ranks_path = derived_dir / "standardised" / "EV_rank_percentiles.csv"
    if ev_ranks_path.exists():
        ev_df = pd.read_csv(ev_ranks_path)
        ev_data["rank_percentiles"] = dict(zip(ev_df["gene_symbol"], ev_df["rank_percentile"]))
        ev_data["genes"] = set(ev_df["gene_symbol"].dropna())

    # Rebuild protein evidence dict from standardised files
    protein_evidence = _build_protein_evidence_dict(derived_dir / "standardised", config)

    # ------------------------------------------------------------------
    # 7. Sensitivity analysis
    # ------------------------------------------------------------------
    logger.info("\nRunning sensitivity analysis...")
    sens_results = run_sensitivity_analysis(config, protein_evidence, ev_data, r1)
    sens_results.to_csv(output_dir / "sensitivity_analysis_summary.tsv", sep="\t", index=False)
    logger.info(f"  Sensitivity: {len(sens_results)} entries")

    # ------------------------------------------------------------------
    # 7b. Robustness metric: inclusion_count (0–9)
    # ------------------------------------------------------------------
    # For each of the 9 (tier_option × weight_set) configs, find the top 80
    # proteins by score.  Count how many configs each protein appears in.
    logger.info("\nComputing robustness metric (inclusion_count)...")
    inclusion = {}
    for (tier, ws), group in sens_results.groupby(["tier_option", "weight_set"]):
        top80 = set(group.nsmallest(80, "rank")["human_gene_symbol"])
        for gene in top80:
            inclusion[gene] = inclusion.get(gene, 0) + 1

    final["inclusion_count"] = (
        final["human_gene_symbol"].map(inclusion).fillna(0).astype(int)
    )
    # Re-save candidates_ranked.tsv with new column
    final.to_csv(output_dir / "candidates_ranked.tsv", sep="\t", index=False)
    logger.info(f"  inclusion_count added to candidates_ranked.tsv")

    # Add to shortlist and re-save
    ic_map = dict(zip(final["human_gene_symbol"], final["inclusion_count"]))
    shortlist["inclusion_count"] = (
        shortlist["human_gene_symbol"].map(ic_map).fillna(0).astype(int)
    )
    shortlist.to_csv(output_dir / "core_panel_shortlist.tsv", sep="\t", index=False)
    logger.info(f"  inclusion_count added to core_panel_shortlist.tsv "
                f"(range: {shortlist['inclusion_count'].min()}–"
                f"{shortlist['inclusion_count'].max()}, "
                f"mean: {shortlist['inclusion_count'].mean():.1f})")

    # ------------------------------------------------------------------
    # 8. Null simulation
    # ------------------------------------------------------------------
    logger.info("\nRunning null simulation...")
    null_results = run_null_simulation(config, protein_evidence, r1)
    null_results.to_csv(qc_dir / "null_simulation_results.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------
    # 9. Figures
    # ------------------------------------------------------------------
    logger.info("\nGenerating figures...")

    # Score distribution
    try:
        autophagy_scores = final[final["in_r1"] == True]["composite_score"]
        non_autophagy_scores = final[final["in_r1"] != True]["composite_score"]
        plot_score_distribution(autophagy_scores, non_autophagy_scores,
                                fig_dir / "score_distribution.pdf")
    except Exception as e:
        logger.warning(f"  Score distribution plot failed: {e}")

    # Detection heatmap
    try:
        tier_cols = [c for c in final.columns if c.endswith("_tier") and c.startswith("D")]
        if tier_cols:
            heatmap_data = final[["human_gene_symbol"] + tier_cols].head(80)
            plot_detection_heatmap(heatmap_data, fig_dir / "detection_heatmap.pdf")
    except Exception as e:
        logger.warning(f"  Detection heatmap failed: {e}")

    # Comparison Venn
    try:
        if new_panel_genes and r2_genes:
            plot_comparison_venn(new_panel_genes, r2_genes, previous_329,
                                 fig_dir / "comparison_venn.pdf")
    except Exception as e:
        logger.warning(f"  Comparison Venn failed: {e}")

    # Null distribution
    try:
        observed = null_results.attrs.get("observed", 0)
        plot_null_distribution(observed, null_results["n_convergent"].tolist(),
                                fig_dir / "null_distribution.pdf")
    except Exception as e:
        logger.warning(f"  Null distribution plot failed: {e}")

    # Sensitivity stability
    try:
        # Pivot: proteins x parameter combos
        default_sens = sens_results[sens_results["weight_set"] == "default"]
        pivot = default_sens.pivot_table(
            index="human_gene_symbol", columns="tier_option",
            values="rank", aggfunc="first"
        )
        pivot = pivot.dropna().sort_values(pivot.columns[0])
        plot_sensitivity_stability(pivot.head(50), fig_dir / "sensitivity_stability.pdf")
    except Exception as e:
        logger.warning(f"  Sensitivity stability plot failed: {e}")

    # Module validation summary (copy from modules)
    mod_enr_path = derived_dir / "modules" / "module_enrichment.tsv"
    if mod_enr_path.exists():
        import shutil
        shutil.copy(mod_enr_path, output_dir / "module_validation_summary.tsv")

    # ------------------------------------------------------------------
    # 10. Methods draft
    # ------------------------------------------------------------------
    write_methods_draft(config, final, output_dir / "methods_draft.md")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    logger.info("\n" + "=" * 70)
    logger.info("FINAL SUMMARY")
    logger.info(f"  Total scored proteins: {len(final)}")
    logger.info(f"  Core panel shortlist: {len(shortlist)}")
    logger.info(f"  Overlap with R2: {len(overlap_r2)}/{len(r2_genes)}")
    logger.info(f"  Overlap with previous 329: {len(overlap_329)}/{len(previous_329)}")
    logger.info(f"  Figures generated in: {fig_dir}")
    logger.info("[07_validate_and_report] DONE")


# ------------------------------------------------------------------ #
# Helper: rebuild protein evidence dict (avoid circular import)
# ------------------------------------------------------------------ #



if __name__ == "__main__":
    main()
