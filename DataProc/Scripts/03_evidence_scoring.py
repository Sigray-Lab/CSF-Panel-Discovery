#!/usr/bin/env python3
"""Step 3: Evidence scoring.

Computes per-protein continuous evidence scores combining:
- Mouse CSF evidence (w1=0.25)
- Human CSF evidence (w2=0.30)
- EV evidence (w3=0.05)
- Brain plausibility (w4=0.10)
- Autophagy membership (w5=0.25)
- Penalties for ambiguity, orthology, plasma contamination

Outputs:
- DerivedData/evidence_scores/evidence_scores.tsv
- DerivedData/evidence_scores/provenance_matrix.tsv
- Outputs/figures/upset_evidence.pdf
"""

import json
import logging
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS_DIR))

from utils.parsers import load_config, resolve_path, parse_r1_reference
from utils.scoring import (
    compute_mouse_csf_evidence, compute_human_csf_evidence,
    compute_ev_evidence, compute_brain_plausibility,
    compute_autophagy_membership, compute_penalties,
    compute_composite_score, categorize_r1_category,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("03_evidence_scoring")


def build_protein_evidence_dict(std_dir: Path, config: dict) -> dict:
    """Merge all standardised datasets into per-protein evidence dictionary.

    Key: human_gene_symbol (uppercase for matching)
    Value: dict with per-dataset detection info
    """
    protein_evidence = {}
    all_datasets = list(config["datasets"].keys())

    for ds_id in all_datasets:
        if ds_id == "EV":
            continue  # handled separately

        ds_path = std_dir / f"{ds_id}_standardised.csv"
        if not ds_path.exists():
            logger.warning(f"  {ds_id}: standardised file not found")
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

            # Per-dataset info
            detected = (row.get("n_detected", 0) > 0
                        or row.get("detectability_tier", "absent") != "absent")
            protein_evidence[gene][ds_id] = {
                "detected": detected,
                "n_detected": int(row.get("n_detected", 0)),
                "pct_detected": float(row.get("pct_detected", 0.0)),
                "median_intensity": float(row.get("median_intensity", 0.0)),
                "detectability_tier": row.get("detectability_tier", "absent"),
            }

            # Aggregate flags (True if True in ANY dataset)
            if row.get("is_ambiguous_group", False):
                protein_evidence[gene]["is_ambiguous_group"] = True
            if row.get("orthology_ambiguous", False):
                protein_evidence[gene]["orthology_ambiguous"] = True
            if row.get("likely_plasma_derived", False):
                protein_evidence[gene]["likely_plasma_derived"] = True

    logger.info(f"Built evidence dict: {len(protein_evidence)} unique proteins")
    return protein_evidence


def build_provenance_matrix(protein_evidence: dict, dataset_ids: list[str]) -> pd.DataFrame:
    """Build protein x dataset detection/tier matrix."""
    records = []
    for gene, evidence in protein_evidence.items():
        row = {"human_gene_symbol": gene}
        for ds_id in dataset_ids:
            ds_info = evidence.get(ds_id, {})
            row[f"{ds_id}_detected"] = ds_info.get("detected", False)
            row[f"{ds_id}_tier"] = ds_info.get("detectability_tier", "absent")
            row[f"{ds_id}_pct_detected"] = ds_info.get("pct_detected", 0.0)
        records.append(row)

    return pd.DataFrame(records).sort_values("human_gene_symbol")


def main():
    config_path = SCRIPTS_DIR / "config.yaml"
    config = load_config(config_path)
    raw_dir = (SCRIPTS_DIR / config["paths"]["raw_dir"]).resolve()
    derived_dir = (SCRIPTS_DIR / config["paths"]["derived_dir"]).resolve()
    output_dir = (SCRIPTS_DIR / config["paths"]["output_dir"]).resolve()
    log_dir = (SCRIPTS_DIR / config["paths"]["log_dir"]).resolve()
    std_dir = derived_dir / "standardised"
    scores_dir = derived_dir / "evidence_scores"
    scores_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    file_handler = logging.FileHandler(log_dir / f"03_evidence_scoring_{timestamp}.log")
    file_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logging.getLogger().addHandler(file_handler)

    scoring_cfg = config["scoring"]
    weights = scoring_cfg["weights"]

    logger.info("=" * 70)
    logger.info("Step 3: Evidence Scoring")
    logger.info(f"Weights: {weights}")
    logger.info("=" * 70)

    # ------------------------------------------------------------------
    # 1. Build protein-level evidence dictionary
    # ------------------------------------------------------------------
    protein_evidence = build_protein_evidence_dict(std_dir, config)

    # ------------------------------------------------------------------
    # 2. Load EV data
    # ------------------------------------------------------------------
    ev_ranks_path = std_dir / "EV_rank_percentiles.csv"
    ev_data = {"rank_percentiles": {}, "genes": set()}
    if ev_ranks_path.exists():
        ev_df = pd.read_csv(ev_ranks_path)
        ev_data["rank_percentiles"] = dict(zip(ev_df["gene_symbol"], ev_df["rank_percentile"]))
        ev_data["genes"] = set(ev_df["gene_symbol"].dropna().astype(str))
        logger.info(f"Loaded EV data: {len(ev_data['genes'])} genes")
    else:
        logger.warning("EV rank percentiles not found")

    # ------------------------------------------------------------------
    # 3. Load R1 reference
    # ------------------------------------------------------------------
    r1_ref = config["references"]["R1"]
    r1_path = resolve_path(raw_dir, r1_ref["file"])
    r1 = parse_r1_reference(r1_path, r1_ref)
    category_modifiers = scoring_cfg["autophagy_category_modifiers"]

    # ------------------------------------------------------------------
    # 4. Score every protein
    # ------------------------------------------------------------------
    logger.info("Computing evidence scores...")
    mouse_csf_ds = scoring_cfg["mouse_csf_datasets"]
    tier_scores = scoring_cfg["human_csf_tier_scores"]
    d10_bonus = scoring_cfg["d10_bonus"]
    penalty_cfg = scoring_cfg["penalties"]

    results = []
    for gene, evidence in protein_evidence.items():
        # Components
        mouse_csf = compute_mouse_csf_evidence(evidence, mouse_csf_ds)
        human_csf = compute_human_csf_evidence(evidence, tier_scores, d10_bonus)
        ev = compute_ev_evidence(gene, ev_data)
        brain = compute_brain_plausibility(evidence)
        autophagy = compute_autophagy_membership(
            gene, r1["genes"], r1["category_map"], category_modifiers
        )

        # Penalties
        penalties = compute_penalties(
            evidence.get("is_ambiguous_group", False),
            evidence.get("orthology_ambiguous", False),
            evidence.get("likely_plasma_derived", False),
            penalty_cfg,
        )

        # Composite
        components = {
            "mouse_csf": mouse_csf,
            "human_csf": human_csf,
            "ev": ev,
            "brain_plausibility": brain,
            "autophagy_membership": autophagy,
        }
        composite = compute_composite_score(components, weights, penalties)

        # Count mouse CSF datasets detected
        n_mouse_ds = sum(
            1 for ds_id in mouse_csf_ds
            if evidence.get(ds_id, {}).get("detected", False)
        )

        # D11/D12/D10/D6/D9 flags
        d11_tier = evidence.get("D11", {}).get("detectability_tier", "absent")
        d12_detected = evidence.get("D12", {}).get("detected", False)
        d10_detected = evidence.get("D10", {}).get("detected", False)
        d6_detected = evidence.get("D6", {}).get("detected", False)
        d9_detected = evidence.get("D9", {}).get("detected", False)

        # Best mouse CSF tier
        mouse_tiers = []
        for ds_id in mouse_csf_ds:
            ds_info = evidence.get(ds_id, {})
            if ds_info.get("detected", False):
                mouse_tiers.append(ds_info.get("detectability_tier", "C"))
        tier_order = {"A": 0, "B": 1, "C": 2, "absent": 3}
        best_mouse_tier = min(mouse_tiers, key=lambda t: tier_order.get(t, 3)) if mouse_tiers else "absent"

        # R1 category
        raw_cat = r1["category_map"].get(gene, "")
        r1_tier = categorize_r1_category(raw_cat)

        results.append({
            "human_gene_symbol": gene,
            "composite_score": round(composite, 6),
            "score_mouse_csf": round(mouse_csf, 6),
            "score_human_csf": round(human_csf, 6),
            "score_ev": round(ev, 6),
            "score_brain": round(brain, 6),
            "score_autophagy": round(autophagy, 6),
            "penalties": round(penalties, 6),
            "n_mouse_csf_datasets": n_mouse_ds,
            "mouse_csf_best_tier": best_mouse_tier,
            "d11_tier": d11_tier,
            "d12_validated": d12_detected,
            "d10_detected": d10_detected,
            "d6_brain_detected": d6_detected,
            "d9_isf_detected": d9_detected,
            "ev_present": gene.upper() in {g.upper() for g in ev_data["genes"]},
            "in_r1": gene.upper() in {g.upper() for g in r1["genes"]},
            "r1_category": raw_cat if raw_cat else None,
            "r1_tier": r1_tier if gene.upper() in {g.upper() for g in r1["genes"]} else None,
            "is_ambiguous_group": evidence.get("is_ambiguous_group", False),
            "orthology_ambiguous": evidence.get("orthology_ambiguous", False),
            "likely_plasma_derived": evidence.get("likely_plasma_derived", False),
        })

    scores_df = pd.DataFrame(results).sort_values("composite_score", ascending=False)
    scores_df = scores_df.reset_index(drop=True)
    scores_df.index.name = "rank"

    logger.info(f"Scored {len(scores_df)} proteins")
    logger.info(f"Score range: {scores_df['composite_score'].min():.4f} - "
                f"{scores_df['composite_score'].max():.4f}")
    logger.info(f"In R1: {scores_df['in_r1'].sum()}")
    logger.info(f"Plasma flagged: {scores_df['likely_plasma_derived'].sum()}")

    # ------------------------------------------------------------------
    # 5. Build provenance matrix
    # ------------------------------------------------------------------
    proteomics_ds = [k for k in config["datasets"].keys() if k != "EV"]
    provenance = build_provenance_matrix(protein_evidence, proteomics_ds)

    # ------------------------------------------------------------------
    # 6. Generate UpSet plot
    # ------------------------------------------------------------------
    try:
        from utils.visualization import plot_upset

        mouse_csf_genes = set(
            scores_df[scores_df["n_mouse_csf_datasets"] > 0]["human_gene_symbol"]
        )
        human_csf_genes = set(
            scores_df[scores_df["d11_tier"] != "absent"]["human_gene_symbol"]
        )
        ev_genes = set(scores_df[scores_df["ev_present"]]["human_gene_symbol"])
        brain_genes = set(scores_df[scores_df["d6_brain_detected"]]["human_gene_symbol"])

        sets = {
            "Mouse CSF": mouse_csf_genes,
            "Human CSF": human_csf_genes,
            "EV": ev_genes,
            "Brain": brain_genes,
            "R1 Autophagy": r1["genes"],
        }

        fig_dir = output_dir / "figures"
        fig_dir.mkdir(parents=True, exist_ok=True)
        plot_upset(sets, "Multi-evidence Overlap", fig_dir / "upset_evidence.pdf")
        logger.info("Generated UpSet plot")
    except Exception as e:
        logger.warning(f"Could not generate UpSet plot: {e}")

    # ------------------------------------------------------------------
    # 7. Save outputs
    # ------------------------------------------------------------------
    scores_path = scores_dir / "evidence_scores.tsv"
    scores_df.to_csv(scores_path, sep="\t", index=True)
    logger.info(f"Evidence scores: {scores_path}")

    provenance_path = scores_dir / "provenance_matrix.tsv"
    provenance.to_csv(provenance_path, sep="\t", index=False)
    logger.info(f"Provenance matrix: {provenance_path}")

    # Quick distribution summary
    logger.info("\nScore distribution (top 20):")
    for _, row in scores_df.head(20).iterrows():
        logger.info(f"  {row['human_gene_symbol']:15s} "
                    f"score={row['composite_score']:.4f} "
                    f"R1={'Y' if row['in_r1'] else 'N'} "
                    f"D11={row['d11_tier']} "
                    f"mouse={row['n_mouse_csf_datasets']}ds")

    logger.info("[03_evidence_scoring] DONE")


if __name__ == "__main__":
    main()
