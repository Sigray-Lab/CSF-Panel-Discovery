#!/usr/bin/env python3
"""Step 1: Extract, QC, and standardise all proteomics datasets.

Per dataset:
1. Parse identifiers (MaxQuant or DIA-NN format)
2. Remove contaminants and flag plasma proteins
3. Apply sample-level QC
4. Compute detection stats and assign detectability tiers (A/B/C)
5. Apply protein-group ambiguity rules
6. Write standardised CSV to DerivedData/standardised/

Outputs:
- DerivedData/standardised/{dataset_id}_standardised.csv (one per dataset)
- DerivedData/standardised/EV_standardised.csv
- QC/sample_missingness_flags.tsv
- QC/contaminant_report.tsv
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

from utils.parsers import (
    load_config, resolve_path, identify_sample_columns, standardize_columns,
    parse_maxquant_fasta, parse_maxquant_preparsed,
    parse_diann_tsv, parse_diann_excel, parse_diann_astral,
    parse_ev_dataset, parse_r1_reference, parse_r2_reference,
)
from utils.qc import (
    remove_contaminants, flag_plasma_proteins, compute_sample_qc,
    detect_empty_columns, compute_intensity_floor, compute_detection_stats,
    assign_tiers_bulk,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    handlers=[logging.StreamHandler()],
)
logger = logging.getLogger("01_extract_and_qc")


# ------------------------------------------------------------------ #
# Parser dispatch
# ------------------------------------------------------------------ #

def parse_dataset(dataset_id: str, ds_config: dict, raw_dir: Path) -> pd.DataFrame:
    """Dispatch to the correct parser based on format type."""
    fmt = ds_config["format"]
    filepath = resolve_path(raw_dir, ds_config["file"])

    if fmt == "maxquant_fasta":
        return parse_maxquant_fasta(filepath, ds_config["sheet"], ds_config)
    elif fmt == "maxquant_preparsed":
        return parse_maxquant_preparsed(filepath, ds_config["sheet"], ds_config)
    elif fmt == "diann_tsv":
        return parse_diann_tsv(filepath, ds_config)
    elif fmt == "diann_excel":
        return parse_diann_excel(filepath, ds_config)
    elif fmt == "diann_astral":
        return parse_diann_astral(filepath, ds_config)
    else:
        raise ValueError(f"Unknown format for {dataset_id}: {fmt}")


# ------------------------------------------------------------------ #
# Main pipeline
# ------------------------------------------------------------------ #

def process_single_dataset(dataset_id: str, ds_config: dict, config: dict,
                            raw_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Process a single proteomics dataset through the full QC pipeline.

    Returns: (standardised_df, sample_qc_df, contaminant_info)
    """
    logger.info(f"{'='*60}")
    logger.info(f"Processing {dataset_id} ({ds_config['species']} {ds_config['biofluid']})")
    logger.info(f"{'='*60}")

    # 1. Parse
    df = parse_dataset(dataset_id, ds_config, raw_dir)
    initial_rows = len(df)
    logger.info(f"  Parsed: {initial_rows} protein groups")

    # 2. Identify sample columns
    sample_cols = identify_sample_columns(df, ds_config)
    logger.info(f"  Sample columns: {len(sample_cols)}")

    # Validate expected sample count
    expected = ds_config.get("n_samples_expected")
    if expected and abs(len(sample_cols) - expected) > expected * 0.20:
        logger.warning(
            f"  WARNING: Expected ~{expected} samples, found {len(sample_cols)}"
        )

    # 3. Remove contaminants (MaxQuant formats)
    contaminant_info = {"dataset_id": dataset_id, "n_removed": 0, "n_before": initial_rows}
    if ds_config["format"].startswith("maxquant"):
        prefixes = config.get("contaminants", {}).get("maxquant_prefixes", ["REV__", "CON__"])
        # Check in the original protein IDs column
        pid_col = ds_config.get("id_columns", {}).get("protein_ids", "group_id")
        check_col = pid_col if pid_col in df.columns else "group_id"
        df, n_removed = remove_contaminants(df, id_col=check_col, prefixes=prefixes)
        contaminant_info["n_removed"] = n_removed
        logger.info(f"  After contaminant removal: {len(df)} rows")

    # 4. Detect and exclude empty sample columns (esp. D9)
    empty_cols = detect_empty_columns(df, sample_cols)
    if empty_cols:
        sample_cols = [c for c in sample_cols if c not in empty_cols]
        logger.info(f"  Excluded {len(empty_cols)} empty columns, {len(sample_cols)} remain")

    # 5. Sample-level QC
    sample_qc = compute_sample_qc(df, sample_cols, dataset_id)
    bad_samples = sample_qc[sample_qc["flagged_high_missingness"]]["sample_name"].tolist()
    if bad_samples:
        logger.warning(f"  Excluding {len(bad_samples)} high-missingness samples")
        sample_cols = [c for c in sample_cols if c not in bad_samples]

    # 6. Replace zeros with NaN for MaxQuant datasets (zeros = missing convention)
    # DIA-NN already uses NaN for missing; still ensure numeric dtype
    for col in sample_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    if ds_config["format"].startswith("maxquant"):
        for col in sample_cols:
            zero_mask = df[col] == 0
            if zero_mask.any():
                df.loc[zero_mask, col] = np.nan

    # 7. Compute detection stats
    df = compute_detection_stats(df, sample_cols)

    # 8. Compute intensity floor and assign detectability tiers
    intensity_floor = compute_intensity_floor(df, sample_cols,
                                               config["detectability"]["intensity_floor_percentile"])
    df["above_intensity_floor"] = df["median_intensity"] > intensity_floor
    df["detectability_tier"] = assign_tiers_bulk(df, sample_cols, intensity_floor,
                                                  config["detectability"])

    tier_counts = df["detectability_tier"].value_counts()
    logger.info(f"  Tier distribution: {tier_counts.to_dict()}")

    # 9. Flag plasma proteins
    df["likely_plasma_derived"] = flag_plasma_proteins(df, "gene_symbol", config)

    # 10. Standardise columns
    df = standardize_columns(df, dataset_id, ds_config["species"], sample_cols)

    # 11. Validate: fail if zero rows
    if len(df) == 0:
        raise ValueError(f"{dataset_id}: zero rows after QC pipeline!")

    logger.info(f"  Final: {len(df)} protein groups, {len(sample_cols)} samples")
    return df, sample_qc, contaminant_info


def main():
    config_path = SCRIPTS_DIR / "config.yaml"
    config = load_config(config_path)
    raw_dir = (SCRIPTS_DIR / config["paths"]["raw_dir"]).resolve()
    derived_dir = (SCRIPTS_DIR / config["paths"]["derived_dir"]).resolve()
    qc_dir = (SCRIPTS_DIR / config["paths"]["qc_dir"]).resolve()
    log_dir = (SCRIPTS_DIR / config["paths"]["log_dir"]).resolve()

    # Set up file logging
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    file_handler = logging.FileHandler(log_dir / f"01_extract_and_qc_{timestamp}.log")
    file_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logging.getLogger().addHandler(file_handler)

    logger.info("=" * 70)
    logger.info("Step 1: Extract, QC, and Standardise")
    logger.info("=" * 70)

    std_dir = derived_dir / "standardised"
    std_dir.mkdir(parents=True, exist_ok=True)

    all_sample_qc = []
    all_contaminant_info = []
    summary = {}

    # ------------------------------------------------------------------ #
    # Process proteomics datasets (D1-D12)
    # ------------------------------------------------------------------ #
    for dataset_id, ds_config in config["datasets"].items():
        if ds_config["format"] == "ev_rank":
            continue  # Handle EV separately

        try:
            df, sample_qc, contam_info = process_single_dataset(
                dataset_id, ds_config, config, raw_dir
            )

            # Save standardised CSV
            output_path = std_dir / f"{dataset_id}_standardised.csv"
            df.to_csv(output_path, index=False)
            logger.info(f"  Saved: {output_path}")

            all_sample_qc.append(sample_qc)
            all_contaminant_info.append(contam_info)
            summary[dataset_id] = {
                "n_proteins": len(df),
                "n_samples": int(df["n_detected"].max()) if "n_detected" in df.columns else 0,
                "tier_distribution": df["detectability_tier"].value_counts().to_dict(),
                "n_plasma_flagged": int(df["likely_plasma_derived"].astype(bool).sum()),
                "n_ambiguous_groups": int(df["is_ambiguous_group"].astype(bool).sum()),
            }

        except Exception as e:
            logger.error(f"FAILED processing {dataset_id}: {e}", exc_info=True)
            raise

    # ------------------------------------------------------------------ #
    # Process EV dataset
    # ------------------------------------------------------------------ #
    ev_config = config["datasets"].get("EV")
    if ev_config:
        logger.info("=" * 60)
        logger.info("Processing EV dataset")
        logger.info("=" * 60)
        ev_path = resolve_path(raw_dir, ev_config["file"])
        ev_data = parse_ev_dataset(ev_path, ev_config)

        # Save as JSON (dict structure) and CSV (precomputed table)
        ev_output = std_dir / "EV_standardised.csv"
        ev_data["precomputed"].to_csv(ev_output, index=False)
        logger.info(f"  Saved EV precomputed: {ev_output}")

        # Save rank percentiles as TSV
        rank_df = pd.DataFrame([
            {"gene_symbol": gene, "rank_percentile": pct}
            for gene, pct in ev_data["rank_percentiles"].items()
        ]).sort_values("rank_percentile")
        rank_output = std_dir / "EV_rank_percentiles.csv"
        rank_df.to_csv(rank_output, index=False)
        logger.info(f"  Saved EV ranks: {rank_output} ({len(rank_df)} genes)")

        summary["EV"] = {
            "n_wt_genes": len(ev_data["raw_wt"]),
            "n_ko_genes": len(ev_data["raw_ko"]),
            "n_unique_genes": len(ev_data["genes"]),
            "n_precomputed": len(ev_data["precomputed"]),
        }

    # ------------------------------------------------------------------ #
    # Write QC outputs
    # ------------------------------------------------------------------ #

    # Sample missingness report
    if all_sample_qc:
        qc_combined = pd.concat(all_sample_qc, ignore_index=True)
        qc_path = qc_dir / "sample_missingness_flags.tsv"
        qc_combined.to_csv(qc_path, sep="\t", index=False)
        logger.info(f"\nSample QC report: {qc_path} ({len(qc_combined)} samples across all datasets)")

    # Contaminant report
    if all_contaminant_info:
        contam_df = pd.DataFrame(all_contaminant_info)
        contam_path = qc_dir / "contaminant_report.tsv"
        contam_df.to_csv(contam_path, sep="\t", index=False)
        logger.info(f"Contaminant report: {contam_path}")

    # ------------------------------------------------------------------ #
    # Summary
    # ------------------------------------------------------------------ #
    logger.info("\n" + "=" * 70)
    logger.info("SUMMARY")
    logger.info("=" * 70)
    for ds_id, info in summary.items():
        logger.info(f"  {ds_id}: {info}")

    # Save summary to log
    summary_path = log_dir / f"01_summary_{timestamp}.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2, default=str)
    logger.info(f"\nSummary saved: {summary_path}")
    logger.info("[01_extract_and_qc] DONE")


if __name__ == "__main__":
    main()
