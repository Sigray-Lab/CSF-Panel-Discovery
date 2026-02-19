"""Quality control utilities: contaminant removal, plasma flagging,
sample QC, intensity floor, and detectability tier assignment."""

import re
import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ------------------------------------------------------------------ #
# Contaminant removal
# ------------------------------------------------------------------ #

def remove_contaminants(df: pd.DataFrame, id_col: str = "group_id",
                        prefixes: list[str] | None = None) -> tuple[pd.DataFrame, int]:
    """Remove MaxQuant decoys (REV__) and common contaminants (CON__).

    Parameters
    ----------
    df : DataFrame with protein data
    id_col : column containing protein identifiers to check
    prefixes : list of prefix strings to remove (default: REV__, CON__)

    Returns
    -------
    (filtered_df, n_removed)
    """
    if prefixes is None:
        prefixes = ["REV__", "CON__"]

    mask = pd.Series(False, index=df.index)
    for prefix in prefixes:
        mask |= df[id_col].astype(str).str.contains(prefix, na=False)

    n_removed = mask.sum()
    if n_removed > 0:
        logger.info(f"  Removed {n_removed} contaminant/decoy entries")

    return df[~mask].copy(), n_removed


# ------------------------------------------------------------------ #
# Plasma protein flagging
# ------------------------------------------------------------------ #

def flag_plasma_proteins(df: pd.DataFrame, gene_col: str = "gene_symbol",
                         config: dict | None = None) -> pd.Series:
    """Return boolean Series flagging likely plasma-derived proteins.

    Matching logic:
    - Exact match against built-in list (case-insensitive)
    - Regex matching for patterns (IGH*, IGL*, IGK*, KRT*)
    - Optional external reference list
    """
    if config is None:
        config = {}

    plasma_cfg = config.get("plasma_proteins", {})
    built_in = set(g.upper() for g in plasma_cfg.get("built_in", []))
    patterns = plasma_cfg.get("patterns", [])
    external_path = plasma_cfg.get("external_list")

    genes = df[gene_col].astype(str).str.strip().str.upper()

    # Exact match
    is_plasma = genes.isin(built_in)

    # Pattern match (regex)
    for pattern in patterns:
        is_plasma |= genes.str.match(pattern, case=False, na=False)

    # External list
    if external_path:
        try:
            ext_df = pd.read_csv(external_path, sep="\t")
            ext_genes = set(ext_df.iloc[:, 0].astype(str).str.upper())
            is_plasma |= genes.isin(ext_genes)
            logger.info(f"  Applied external plasma list: {len(ext_genes)} genes")
        except Exception as e:
            logger.warning(f"  Could not load external plasma list: {e}")

    n_flagged = is_plasma.sum()
    logger.info(f"  Flagged {n_flagged} likely plasma-derived proteins")
    return is_plasma


# ------------------------------------------------------------------ #
# Sample-level QC
# ------------------------------------------------------------------ #

def compute_sample_qc(df: pd.DataFrame, sample_cols: list[str],
                       dataset_id: str) -> pd.DataFrame:
    """Compute per-sample missingness statistics.

    Returns DataFrame with columns:
    dataset_id, sample_name, n_proteins_total, n_detected, pct_missing,
    flagged_high_missingness
    """
    n_total = len(df)
    records = []
    for col in sample_cols:
        n_detected = df[col].notna().sum()
        pct_missing = 1.0 - (n_detected / max(n_total, 1))
        records.append({
            "dataset_id": dataset_id,
            "sample_name": col,
            "n_proteins_total": n_total,
            "n_detected": n_detected,
            "pct_missing": round(pct_missing, 4),
            "flagged_high_missingness": pct_missing > 0.80,
        })

    qc_df = pd.DataFrame(records)
    n_flagged = qc_df["flagged_high_missingness"].sum()
    if n_flagged > 0:
        logger.warning(f"  {dataset_id}: {n_flagged}/{len(sample_cols)} samples flagged "
                       f"with >80% missingness")
    return qc_df


def detect_empty_columns(df: pd.DataFrame, sample_cols: list[str],
                          threshold: float = 0.99) -> list[str]:
    """Identify sample columns that are >threshold fraction empty.

    Used for D9 (2 completely empty columns).
    """
    empty_cols = []
    n_total = len(df)
    for col in sample_cols:
        pct_missing = df[col].isna().sum() / max(n_total, 1)
        if pct_missing >= threshold:
            empty_cols.append(col)

    if empty_cols:
        logger.info(f"  Detected {len(empty_cols)} empty columns (>={threshold*100}% missing)")
    return empty_cols


# ------------------------------------------------------------------ #
# Detection statistics
# ------------------------------------------------------------------ #

def compute_intensity_floor(df: pd.DataFrame, sample_cols: list[str],
                              percentile: int = 10) -> float:
    """Compute dataset-specific intensity floor.

    Returns the given percentile of all non-NaN intensity values
    across all proteins and samples.
    """
    all_values = df[sample_cols].values.flatten()
    all_values = all_values[~np.isnan(all_values)]
    all_values = all_values[all_values > 0]

    if len(all_values) == 0:
        logger.warning("  No non-zero intensity values found; floor = 0")
        return 0.0

    floor = float(np.percentile(all_values, percentile))
    logger.info(f"  Intensity floor ({percentile}th percentile): {floor:.4g}")
    return floor


def compute_detection_stats(df: pd.DataFrame, sample_cols: list[str]) -> pd.DataFrame:
    """Compute per-protein detection statistics.

    Adds columns: n_detected, pct_detected, median_intensity.
    """
    df = df.copy()
    # Ensure all sample columns are numeric
    for col in sample_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    intensity_matrix = df[sample_cols]
    df["n_detected"] = intensity_matrix.notna().sum(axis=1).astype(int)
    df["pct_detected"] = df["n_detected"] / max(len(sample_cols), 1)
    df["median_intensity"] = intensity_matrix.median(axis=1, skipna=True)
    return df


def assign_detectability_tier(n_detected: int, pct_detected: float,
                                median_intensity: float,
                                intensity_floor: float,
                                n_samples: int,
                                detect_config: dict) -> str:
    """Assign detectability tier (A/B/C) per protein per dataset.

    Uses small_dataset_threshold to decide between small/large rules.
    """
    threshold = detect_config.get("small_dataset_threshold", 20)
    is_small = n_samples <= threshold
    above_floor = (not np.isnan(median_intensity)) and (median_intensity > intensity_floor)
    tiers = detect_config.get("tiers", {})

    # Check tier A
    tier_a = tiers.get("A", {})
    if is_small:
        a_rules = tier_a.get("small", {})
        if n_detected >= a_rules.get("min_detected", 2) and (not a_rules.get("above_floor", True) or above_floor):
            return "A"
    else:
        a_rules = tier_a.get("large", {})
        if pct_detected >= a_rules.get("min_fraction", 0.10) and (not a_rules.get("above_floor", True) or above_floor):
            return "A"

    # Check tier B
    tier_b = tiers.get("B", {})
    if is_small:
        b_rules = tier_b.get("small", {})
        if n_detected >= b_rules.get("min_detected", 1) and (not b_rules.get("above_floor", True) or above_floor):
            return "B"
    else:
        b_rules = tier_b.get("large", {})
        meets_fraction = pct_detected >= b_rules.get("min_fraction", 0.02)
        meets_floor = not b_rules.get("above_floor", False) or above_floor
        # Large dataset tier B: fraction OR above_floor (OR logic per spec)
        if meets_fraction or above_floor:
            return "B"

    # Tier C: any detection
    tier_c = tiers.get("C", {})
    if is_small:
        c_rules = tier_c.get("small", {})
        if n_detected >= c_rules.get("min_detected", 1):
            return "C"
    else:
        c_rules = tier_c.get("large", {})
        if n_detected >= c_rules.get("min_detected", 1):
            return "C"

    return "absent"


def assign_tiers_bulk(df: pd.DataFrame, sample_cols: list[str],
                       intensity_floor: float, detect_config: dict) -> pd.Series:
    """Assign detectability tiers for all proteins in a DataFrame."""
    n_samples = len(sample_cols)
    return df.apply(
        lambda row: assign_detectability_tier(
            n_detected=int(row.get("n_detected", 0)),
            pct_detected=float(row.get("pct_detected", 0.0)),
            median_intensity=float(row.get("median_intensity", np.nan)),
            intensity_floor=intensity_floor,
            n_samples=n_samples,
            detect_config=detect_config,
        ),
        axis=1,
    )
