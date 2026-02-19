#!/usr/bin/env python3
"""Step 8: Cross-check pipeline candidates against APP knock-in AD mouse model data.

Parses Table S3 from the AppNL-G-F study (CSF proteomics from APP knock-in mice),
maps mouse gene symbols to human orthologs, and cross-references against the
pipeline's ranked candidate list to identify disease-relevant changes.

Three comparisons are analysed:
  1. AppNL-G-F/NL-G-F vs Appwt/wt   (full knock-in vs wild-type)
  2. AppNL-F/NL-F vs Appwt/wt        (Swedish + Iberian mutations vs WT)
  3. AppNL-F/NL-F vs AppNL-G-F/NL-G-F (partial vs full knock-in)

Outputs:
- Outputs/ad_model_crosscheck.tsv      (annotated overlap table)
- QC/ad_model_crosscheck_summary.txt   (human-readable summary)
"""

import logging
import re
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

SCRIPTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS_DIR))

from utils.parsers import load_config, resolve_path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("08_ad_model_crosscheck")


# ── Parsing ──────────────────────────────────────────────────────────────


def parse_appnlgf_sheet(filepath: Path, sheet_name: str) -> pd.DataFrame:
    """Parse one sheet from the AppNL-G-F Excel file.

    Extracts mouse gene symbols from UniProt FASTA-header entries using
    the GN= tag.  Returns one row per unique gene symbol with ratio and
    p-value.

    Parameters
    ----------
    filepath : Path
        Path to the Excel workbook.
    sheet_name : str
        Sheet name to read.

    Returns
    -------
    pd.DataFrame
        Columns: mouse_gene, ratio, pvalue, avg_wt, avg_mutant, uniprot_entry
    """
    df = pd.read_excel(filepath, sheet_name=sheet_name)
    cols = df.columns.tolist()

    records = []
    for _, row in df.iterrows():
        entry = str(row.iloc[0])
        # Extract all GN= tags (may have multiple from semicolon-separated entries)
        gene_symbols = re.findall(r"GN=(\w+)", entry)
        if not gene_symbols:
            continue
        # Take the first unique gene symbol (multi-accession rows map to same gene)
        mouse_gene = gene_symbols[0]

        # Identify numeric columns by position
        # Ratio is second-to-last numeric column, p-value is last
        ratio = row.iloc[-2]
        pvalue = row.iloc[-1]
        # Averages: the two columns before ratio
        avg_wt = row.iloc[-4]
        avg_mutant = row.iloc[-3]

        records.append({
            "mouse_gene": mouse_gene,
            "ratio": ratio,
            "pvalue": pvalue,
            "avg_wt": avg_wt,
            "avg_mutant": avg_mutant,
            "uniprot_entry": entry[:120],  # truncate for readability
        })

    result = pd.DataFrame(records)
    # Deduplicate: if same gene appears multiple times, keep first occurrence
    result = result.drop_duplicates(subset="mouse_gene", keep="first")
    logger.info(f"  Parsed {len(result)} unique genes from '{sheet_name}'")
    return result


# ── Orthology mapping ────────────────────────────────────────────────────


def load_orthology_map(cache_path: Path) -> dict:
    """Load orthology cache and build mouse→human gene symbol mapping.

    Parameters
    ----------
    cache_path : Path
        Path to orthology_cache.tsv.

    Returns
    -------
    dict
        {mouse_gene_lower: human_gene_symbol}
    """
    cache = pd.read_csv(cache_path, sep="\t")
    mapping = {}
    for _, row in cache.iterrows():
        mouse = str(row["input_symbol"]).lower()
        human = str(row["output_symbol"]).strip()
        if human and human != "nan":
            mapping[mouse] = human
    logger.info(f"  Loaded {len(mapping)} orthology mappings from cache")
    return mapping


def map_mouse_to_human(mouse_genes: list, ortho_map: dict) -> dict:
    """Map mouse gene symbols to human orthologs.

    Strategy:
    1. Look up in orthology cache (case-insensitive)
    2. Fall back to uppercase (many genes share symbol across species)

    Parameters
    ----------
    mouse_genes : list
        Mouse gene symbols from the AppNL-G-F table.
    ortho_map : dict
        {mouse_gene_lower: human_gene_symbol} from orthology cache.

    Returns
    -------
    dict
        {original_mouse_gene: human_gene_symbol}
    """
    result = {}
    n_cache = 0
    n_uppercase = 0
    n_unmapped = 0

    for mg in mouse_genes:
        key = mg.lower()
        if key in ortho_map:
            result[mg] = ortho_map[key]
            n_cache += 1
        else:
            # Fallback: uppercase mouse gene = human gene
            result[mg] = mg.upper()
            n_uppercase += 1

    logger.info(
        f"  Mapped {n_cache} via cache, {n_uppercase} via uppercase fallback"
    )
    return result


# ── Cross-reference ──────────────────────────────────────────────────────


def classify_direction(ratio: float, pvalue: float, threshold: float) -> str:
    """Classify fold-change direction with significance."""
    if pd.isna(ratio) or pd.isna(pvalue):
        return "NA"
    if pvalue >= threshold:
        return "ns"
    return "up" if ratio > 1.0 else "down"


def crosscheck(
    candidates: pd.DataFrame,
    ad_sheets: dict,
    gene_map: dict,
    sig_threshold: float = 0.05,
) -> pd.DataFrame:
    """Cross-reference pipeline candidates against AD model data.

    Parameters
    ----------
    candidates : pd.DataFrame
        Full candidates_ranked.tsv.
    ad_sheets : dict
        {label: DataFrame} from parse_appnlgf_sheet.
    gene_map : dict
        {mouse_gene: human_gene_symbol}.
    sig_threshold : float
        P-value threshold for significance.

    Returns
    -------
    pd.DataFrame
        Annotated overlap table.
    """
    # Build human_gene → AD data mapping for each comparison
    ad_by_human = {}  # {label: {human_gene: {ratio, pvalue, direction}}}

    for label, ad_df in ad_sheets.items():
        hmap = {}
        for _, row in ad_df.iterrows():
            human = gene_map.get(row["mouse_gene"])
            if human:
                hmap[human] = {
                    "ratio": row["ratio"],
                    "pvalue": row["pvalue"],
                    "direction": classify_direction(
                        row["ratio"], row["pvalue"], sig_threshold
                    ),
                    "avg_wt": row["avg_wt"],
                    "avg_mutant": row["avg_mutant"],
                }
        ad_by_human[label] = hmap
        logger.info(f"  {label}: {len(hmap)} proteins mapped to human symbols")

    # Select columns to keep from candidates
    keep_cols = [
        "rank", "human_gene_symbol", "composite_score",
        "n_mouse_csf_datasets", "mouse_csf_best_tier", "d11_tier",
        "in_r1", "r1_category", "inclusion_count",
    ]
    available = [c for c in keep_cols if c in candidates.columns]

    # Find all human genes present in ANY comparison
    all_ad_human = set()
    for hmap in ad_by_human.values():
        all_ad_human.update(hmap.keys())

    # Filter candidates to those present in AD data
    overlap_mask = candidates["human_gene_symbol"].isin(all_ad_human)
    overlap = candidates.loc[overlap_mask, available].copy()

    # Add AD data columns for each comparison
    for label, hmap in ad_by_human.items():
        overlap[f"ratio_{label}"] = overlap["human_gene_symbol"].map(
            lambda g, h=hmap: h.get(g, {}).get("ratio")
        )
        overlap[f"pvalue_{label}"] = overlap["human_gene_symbol"].map(
            lambda g, h=hmap: h.get(g, {}).get("pvalue")
        )
        overlap[f"direction_{label}"] = overlap["human_gene_symbol"].map(
            lambda g, h=hmap: h.get(g, {}).get("direction", "NA")
        )

    # Summary columns
    comparison_labels = list(ad_sheets.keys())
    sig_cols = [f"pvalue_{l}" for l in comparison_labels]
    dir_cols = [f"direction_{l}" for l in comparison_labels]

    overlap["n_comparisons_significant"] = overlap[sig_cols].apply(
        lambda row: sum(
            1 for v in row if pd.notna(v) and v < sig_threshold
        ),
        axis=1,
    )

    def _consistent(row):
        dirs = [
            row[c] for c in dir_cols
            if row.get(c) not in ("ns", "NA", None) and pd.notna(row.get(c))
        ]
        if len(dirs) < 2:
            return True  # trivially consistent (0 or 1 significant)
        return len(set(dirs)) == 1

    overlap["consistent_direction"] = overlap.apply(_consistent, axis=1)

    # Sort by pipeline rank
    overlap = overlap.sort_values("rank").reset_index(drop=True)

    return overlap


# ── Summary report ───────────────────────────────────────────────────────


def write_summary(
    overlap: pd.DataFrame,
    ad_sheets: dict,
    gene_map: dict,
    candidates: pd.DataFrame,
    shortlist: pd.DataFrame,
    output_path: Path,
    sig_threshold: float,
):
    """Write a human-readable summary of the cross-check."""
    n_candidates_total = len(candidates)
    n_core = int((candidates["in_r1"] == True).sum())
    shortlist_genes = set(shortlist["human_gene_symbol"])

    # Count mapped genes per comparison
    all_ad_human = set()
    for hmap_label, ad_df in ad_sheets.items():
        mapped = {gene_map.get(g) for g in ad_df["mouse_gene"] if gene_map.get(g)}
        all_ad_human.update(mapped)

    overlap_core = overlap[overlap["in_r1"] == True]
    overlap_shortlist = overlap[
        overlap["human_gene_symbol"].isin(shortlist_genes)
    ]

    lines = [
        "=" * 70,
        "AD MODEL CROSS-CHECK SUMMARY",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "=" * 70,
        "",
        "SOURCE: Table S3 from APP knock-in mouse study",
        "  AppNL-G-F/NL-G-F = Swedish + Beyreuther/Iberian + Arctic mutations",
        "  AppNL-F/NL-F = Swedish + Beyreuther/Iberian mutations",
        "  Appwt/wt = wild-type controls",
        "",
        "─" * 50,
        "DATASET SIZES",
        "─" * 50,
    ]

    for label, ad_df in ad_sheets.items():
        n_genes = len(ad_df)
        n_mapped = sum(1 for g in ad_df["mouse_gene"] if gene_map.get(g))
        lines.append(f"  {label}: {n_genes} proteins, {n_mapped} mapped to human")

    lines += [
        "",
        "─" * 50,
        "OVERLAP WITH PIPELINE CANDIDATES",
        "─" * 50,
        f"  Total candidates in pipeline:  {n_candidates_total}",
        f"  Core candidates (in_r1=True):  {n_core}",
        f"  Top 80 shortlist:              {len(shortlist_genes)}",
        f"  AD model unique human genes:   {len(all_ad_human)}",
        "",
        f"  Overlap with ALL candidates:   {len(overlap)}",
        f"  Overlap with core (in_r1):     {len(overlap_core)}",
        f"  Overlap with top 80 shortlist: {len(overlap_shortlist)}",
        "",
    ]

    # Significant hits in shortlist
    sig_in_shortlist = overlap_shortlist[
        overlap_shortlist["n_comparisons_significant"] > 0
    ]

    lines += [
        "─" * 50,
        f"SIGNIFICANT HITS IN TOP 80 SHORTLIST (p < {sig_threshold})",
        "─" * 50,
    ]

    if len(sig_in_shortlist) > 0:
        lines.append(
            f"{'Gene':<12} {'Rank':>5} {'Score':>6} "
            f"{'#Sig':>4} {'Consist':>7}  Details"
        )
        lines.append("-" * 75)

        for _, row in sig_in_shortlist.sort_values("rank").iterrows():
            # Build detail string
            details = []
            for label in ad_sheets.keys():
                d = row.get(f"direction_{label}", "NA")
                p = row.get(f"pvalue_{label}")
                r = row.get(f"ratio_{label}")
                if d not in ("NA", None) and pd.notna(d):
                    p_str = f"{p:.4f}" if pd.notna(p) else "NA"
                    r_str = f"{r:.3f}" if pd.notna(r) else "NA"
                    details.append(f"{label[:15]}:{d}(r={r_str},p={p_str})")

            lines.append(
                f"{row['human_gene_symbol']:<12} "
                f"{int(row['rank']):>5} "
                f"{row['composite_score']:>6.3f} "
                f"{int(row['n_comparisons_significant']):>4} "
                f"{'Yes' if row['consistent_direction'] else 'No':>7}  "
                f"{'; '.join(details)}"
            )
    else:
        lines.append("  None found.")

    # All overlapping proteins (even non-significant)
    lines += [
        "",
        "─" * 50,
        "ALL OVERLAPPING PROTEINS (top 80 shortlist, sorted by rank)",
        "─" * 50,
    ]

    if len(overlap_shortlist) > 0:
        lines.append(
            f"{'Gene':<12} {'Rank':>5} {'Score':>6} {'Category':<22} "
            f"{'#Sig':>4}"
        )
        lines.append("-" * 60)
        for _, row in overlap_shortlist.sort_values("rank").iterrows():
            cat = str(row.get("r1_category", ""))[:20]
            lines.append(
                f"{row['human_gene_symbol']:<12} "
                f"{int(row['rank']):>5} "
                f"{row['composite_score']:>6.3f} "
                f"{cat:<22} "
                f"{int(row['n_comparisons_significant']):>4}"
            )

    lines += ["", "=" * 70, "END OF REPORT", "=" * 70, ""]

    output_path.write_text("\n".join(lines))
    logger.info(f"Summary written to {output_path}")


# ── Main ─────────────────────────────────────────────────────────────────


def main():
    logger.info("=" * 60)
    logger.info("STEP 8: AD Model Cross-Check")
    logger.info("=" * 60)

    # ── Load config ──
    config = load_config(SCRIPTS_DIR / "config.yaml")
    raw_dir = resolve_path(SCRIPTS_DIR, config["paths"]["raw_dir"])
    derived_dir = resolve_path(SCRIPTS_DIR, config["paths"]["derived_dir"])
    output_dir = resolve_path(SCRIPTS_DIR, config["paths"]["output_dir"])
    qc_dir = resolve_path(SCRIPTS_DIR, config["paths"]["qc_dir"])

    ad_config = config.get("ad_model")
    if ad_config is None:
        logger.error("No 'ad_model' section in config.yaml — nothing to do.")
        return

    sig_threshold = ad_config.get("significance_threshold", 0.05)

    # ── Parse AppNL-G-F Excel ──
    ad_filepath = raw_dir / ad_config["file"]
    if not ad_filepath.exists():
        logger.error(f"AD model file not found: {ad_filepath}")
        return

    logger.info(f"\nParsing AD model data from: {ad_filepath.name}")
    ad_sheets = {}
    for comp in ad_config["comparisons"]:
        sheet_name = comp["sheet"]
        label = comp["label"]
        ad_sheets[label] = parse_appnlgf_sheet(ad_filepath, sheet_name)

    # ── Load orthology cache ──
    logger.info("\nLoading orthology mappings...")
    cache_path = derived_dir / "orthology_cache.tsv"
    ortho_map = load_orthology_map(cache_path)

    # Collect all mouse genes across sheets
    all_mouse_genes = set()
    for ad_df in ad_sheets.values():
        all_mouse_genes.update(ad_df["mouse_gene"].tolist())

    gene_map = map_mouse_to_human(list(all_mouse_genes), ortho_map)

    # ── Load pipeline candidates ──
    logger.info("\nLoading pipeline candidates...")
    candidates_path = output_dir / "candidates_ranked.tsv"
    candidates = pd.read_csv(candidates_path, sep="\t")
    logger.info(f"  {len(candidates)} total ranked proteins")

    shortlist_path = output_dir / "core_panel_shortlist.tsv"
    shortlist = pd.read_csv(shortlist_path, sep="\t")
    logger.info(f"  {len(shortlist)} shortlisted proteins")

    # ── Cross-reference ──
    logger.info("\nCross-referencing...")
    overlap = crosscheck(candidates, ad_sheets, gene_map, sig_threshold)
    logger.info(f"\n  Total overlapping proteins: {len(overlap)}")

    # Filter view: core candidates only
    overlap_core = overlap[overlap["in_r1"] == True]
    logger.info(f"  Overlapping core candidates (in_r1): {len(overlap_core)}")

    sig_any = overlap[overlap["n_comparisons_significant"] > 0]
    logger.info(f"  Significant in ≥1 comparison: {len(sig_any)}")

    # ── Save outputs ──
    out_path = output_dir / "ad_model_crosscheck.tsv"
    overlap.to_csv(out_path, sep="\t", index=False)
    logger.info(f"\nSaved cross-check table: {out_path}")

    # ── Write summary ──
    summary_path = qc_dir / "ad_model_crosscheck_summary.txt"
    write_summary(
        overlap, ad_sheets, gene_map, candidates, shortlist,
        summary_path, sig_threshold,
    )

    # ── Print key findings ──
    logger.info("\n" + "─" * 50)
    logger.info("KEY FINDINGS")
    logger.info("─" * 50)

    shortlist_genes = set(shortlist["human_gene_symbol"])
    overlap_shortlist = overlap[
        overlap["human_gene_symbol"].isin(shortlist_genes)
    ]
    sig_shortlist = overlap_shortlist[
        overlap_shortlist["n_comparisons_significant"] > 0
    ]

    logger.info(
        f"  {len(overlap_shortlist)} of top 80 shortlist proteins found in AD model data"
    )
    logger.info(
        f"  {len(sig_shortlist)} show significant changes (p < {sig_threshold})"
    )

    if len(sig_shortlist) > 0:
        logger.info("\n  Significant shortlist hits:")
        for _, row in sig_shortlist.sort_values("rank").iterrows():
            n_sig = int(row["n_comparisons_significant"])
            logger.info(
                f"    {row['human_gene_symbol']:<10} rank={int(row['rank']):>3}  "
                f"sig_in={n_sig}/3 comparisons"
            )

    logger.info("\n✓ Step 8 complete.")


if __name__ == "__main__":
    main()
