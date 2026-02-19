#!/usr/bin/env python3
"""Step 2: Orthology mapping (mouse -> human gene symbols).

1. Pre-curate R1 autophagy genes for ambiguous orthology
2. Map all mouse gene symbols to human orthologs via gprofiler
3. Cross-validate ambiguous cases with mygene
4. Cache all mappings
5. Apply mapping to all mouse datasets
6. Human datasets: gene_symbol = human_gene_symbol

Outputs:
- Updated DerivedData/standardised/*.csv with human_gene_symbol column
- DerivedData/orthology_cache.tsv
- QC/R1_orthology_ambiguous_cases.tsv
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
from utils.orthology import (
    map_mouse_to_human_gprofiler, cross_validate_mygene,
    load_orthology_cache, save_orthology_cache,
    resolve_orthology, precurate_r1_orthology,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("02_orthology_mapping")


def apply_orthology_to_dataset(df: pd.DataFrame, cache: pd.DataFrame,
                                policy: str) -> pd.DataFrame:
    """Apply orthology mapping to a mouse dataset.

    For 1:1 mappings: set human_gene_symbol directly.
    For 1:many with expand policy: duplicate rows for each ortholog.
    """
    results = []

    for idx, row in df.iterrows():
        gene = row.get("gene_symbol")
        species = row.get("species", "mouse")

        if pd.isna(gene) or str(gene).strip() == "":
            row_copy = row.copy()
            row_copy["human_gene_symbol"] = np.nan
            row_copy["orthology_ambiguous"] = False
            results.append(row_copy)
            continue

        mappings = resolve_orthology(str(gene), species, cache, policy)

        for mapping in mappings:
            row_copy = row.copy()
            row_copy["human_gene_symbol"] = mapping["human_gene_symbol"]
            row_copy["orthology_ambiguous"] = mapping["orthology_ambiguous"]
            results.append(row_copy)

    result_df = pd.DataFrame(results)
    return result_df.reset_index(drop=True)


def main():
    config_path = SCRIPTS_DIR / "config.yaml"
    config = load_config(config_path)
    raw_dir = (SCRIPTS_DIR / config["paths"]["raw_dir"]).resolve()
    derived_dir = (SCRIPTS_DIR / config["paths"]["derived_dir"]).resolve()
    qc_dir = (SCRIPTS_DIR / config["paths"]["qc_dir"]).resolve()
    log_dir = (SCRIPTS_DIR / config["paths"]["log_dir"]).resolve()
    std_dir = derived_dir / "standardised"

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    file_handler = logging.FileHandler(log_dir / f"02_orthology_{timestamp}.log")
    file_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logging.getLogger().addHandler(file_handler)

    policy = config.get("protein_groups", {}).get("ambiguous_handling", "expand")
    cache_path = derived_dir / "orthology_cache.tsv"

    logger.info("=" * 70)
    logger.info("Step 2: Orthology Mapping")
    logger.info("=" * 70)

    # ------------------------------------------------------------------
    # 1. Load existing cache
    # ------------------------------------------------------------------
    cache = load_orthology_cache(cache_path)

    # ------------------------------------------------------------------
    # 2. Collect all unique mouse gene symbols
    # ------------------------------------------------------------------
    mouse_datasets = ["D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9"]
    mouse_genes = set()

    for ds_id in mouse_datasets:
        ds_path = std_dir / f"{ds_id}_standardised.csv"
        if not ds_path.exists():
            logger.warning(f"  {ds_id} standardised file not found, skipping")
            continue
        df = pd.read_csv(ds_path)
        genes = df["gene_symbol"].dropna().astype(str).str.strip()
        genes = set(g for g in genes if g and g != "nan")
        mouse_genes.update(genes)
        logger.info(f"  {ds_id}: {len(genes)} unique gene symbols")

    logger.info(f"Total unique mouse gene symbols: {len(mouse_genes)}")

    # ------------------------------------------------------------------
    # 3. Pre-curate R1 autophagy genes
    # ------------------------------------------------------------------
    r1_ref = config["references"]["R1"]
    r1_path = resolve_path(raw_dir, r1_ref["file"])
    r1 = parse_r1_reference(r1_path, r1_ref)
    logger.info(f"R1: {len(r1['genes'])} unique autophagy genes")

    # ------------------------------------------------------------------
    # 4. Map mouse genes to human via gprofiler
    # ------------------------------------------------------------------
    already_cached = set()
    if not cache.empty and "input_symbol" in cache.columns:
        already_cached = set(cache["input_symbol"].astype(str).str.upper())

    genes_to_map = [g for g in mouse_genes
                    if g.upper() not in already_cached]

    if genes_to_map:
        logger.info(f"Mapping {len(genes_to_map)} new genes via gprofiler "
                    f"({len(already_cached)} already cached)")
        new_mappings = map_mouse_to_human_gprofiler(genes_to_map)

        # Cross-validate ambiguous cases with mygene
        if not new_mappings.empty:
            ambiguous = new_mappings[new_mappings["orthology_type"] == "1:many"]
            if len(ambiguous) > 0:
                logger.info(f"Cross-validating {len(ambiguous)} ambiguous mappings with mygene")
                pairs = list(zip(ambiguous["input_symbol"], ambiguous["output_symbol"]))
                validation = cross_validate_mygene(pairs)
                if not validation.empty:
                    # Log confirmation rate
                    n_confirmed = validation["mygene_confirmed"].sum()
                    logger.info(f"  mygene confirmed {n_confirmed}/{len(validation)} mappings")

            # Update cache
            cache = pd.concat([cache, new_mappings], ignore_index=True)
            cache = cache.drop_duplicates(subset=["input_symbol", "output_symbol"])
            save_orthology_cache(cache, cache_path)
    else:
        logger.info("All genes already cached, no new API calls needed")

    # ------------------------------------------------------------------
    # 5. R1 pre-curation: flag ambiguous autophagy gene mappings
    # ------------------------------------------------------------------
    # Map R1 human genes to mouse to identify paralog families
    r1_ortho_report = precurate_r1_orthology(r1["genes"], cache)
    if not r1_ortho_report.empty:
        r1_qc_path = qc_dir / "R1_orthology_ambiguous_cases.tsv"
        r1_ortho_report.to_csv(r1_qc_path, sep="\t", index=False)
        logger.info(f"R1 ambiguous cases: {r1_qc_path} ({len(r1_ortho_report)} genes)")

    # ------------------------------------------------------------------
    # 6. Apply orthology to all mouse datasets
    # ------------------------------------------------------------------
    for ds_id in mouse_datasets:
        ds_path = std_dir / f"{ds_id}_standardised.csv"
        if not ds_path.exists():
            continue

        logger.info(f"Applying orthology to {ds_id}...")
        df = pd.read_csv(ds_path)
        n_before = len(df)

        df = apply_orthology_to_dataset(df, cache, policy)
        n_after = len(df)

        n_mapped = df["human_gene_symbol"].notna().sum()
        n_ambiguous = df["orthology_ambiguous"].sum() if "orthology_ambiguous" in df.columns else 0

        logger.info(f"  {ds_id}: {n_before} -> {n_after} rows "
                    f"({n_mapped} mapped, {n_ambiguous} ambiguous)")

        df.to_csv(ds_path, index=False)

    # ------------------------------------------------------------------
    # 7. Human datasets: gene_symbol = human_gene_symbol
    # ------------------------------------------------------------------
    for ds_id in ["D10", "D11", "D12"]:
        ds_path = std_dir / f"{ds_id}_standardised.csv"
        if not ds_path.exists():
            continue

        logger.info(f"Setting human_gene_symbol for {ds_id} (same species)")
        df = pd.read_csv(ds_path)
        df["human_gene_symbol"] = df["gene_symbol"]
        df["orthology_ambiguous"] = False
        df.to_csv(ds_path, index=False)
        logger.info(f"  {ds_id}: {df['human_gene_symbol'].notna().sum()} genes set")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    n_total_cached = len(cache)
    n_1to1 = len(cache[cache.get("orthology_type", pd.Series()) == "1:1"]) if "orthology_type" in cache.columns else 0
    n_1tomany = len(cache[cache.get("orthology_type", pd.Series()) == "1:many"]) if "orthology_type" in cache.columns else 0

    logger.info("\n" + "=" * 70)
    logger.info("SUMMARY")
    logger.info(f"  Cache: {n_total_cached} total mappings")
    logger.info(f"  1:1 mappings: {n_1to1}")
    logger.info(f"  1:many mappings: {n_1tomany}")
    logger.info(f"  Policy: {policy}")
    logger.info("[02_orthology_mapping] DONE")


if __name__ == "__main__":
    main()
