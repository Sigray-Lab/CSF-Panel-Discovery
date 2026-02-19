#!/usr/bin/env python3
"""Step 10: Supplementary gene scoring (R1b).

Scores supplementary curated gene lists against existing pipeline evidence.
For genes already in evidence_scores.tsv, recalculates the autophagy
component (was 0 — not in R1) and recomputes composite scores. Reports
genes not detected in any dataset separately.

Reusable: parameterised via config.yaml supplementary_gene_lists section.
Add new entries there as curated lists grow over time.

Outputs per list (in Outputs/):
  supplementary_{label}_rescored.tsv        — all found genes, old vs new scores
  supplementary_{label}_core_qualified.tsv  — subset passing all 4 hard gates
  supplementary_{label}_not_detected.tsv    — genes absent from all datasets
  supplementary_{label}_combined_ranking.tsv— R1 core + R1b qualified, merged
  supplementary_{label}_report.txt          — human-readable summary
"""

import logging
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS_DIR))

from utils.parsers import load_config, resolve_path
from utils.scoring import compute_composite_score

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("10_supplementary_gene_scoring")


# ------------------------------------------------------------------ #
# HADb category → scoring tier mapping
# ------------------------------------------------------------------ #

def _map_hadb_category(hadb_category: str) -> str:
    """Map HADb Core 'Category' string to a pipeline scoring tier."""
    cat = hadb_category.lower()
    if "mitophagy" in cat:
        return "mitophagy"
    if any(kw in cat for kw in ["fusion", "tethering", "docking"]):
        return "docking"
    if any(kw in cat for kw in ["regulatory", "regulator", "signaling"]):
        return "upstream_regulators"
    if "lysosom" in cat:
        return "lysosomal"
    # CMA, biogenesis, receptor, scaffold, PI3K → core machinery
    return "core_machinery"


# ------------------------------------------------------------------ #
# Database Sources keywords → scoring tier
# ------------------------------------------------------------------ #

_SOURCE_TIER_RULES = [
    # (keyword_in_db_sources, tier, modifier) — checked in order
    ("HADb_core", "core_machinery", 1.0),
    ("GO_lysosome", "lysosomal", 1.0),
    ("KEGG_Lysosome", "lysosomal", 1.0),
    ("Reactome_Lysosome_Biogenesis", "lysosomal", 1.0),
    ("GO_VATPase", "lysosomal", 1.0),
    ("GO_mitophagy", "mitophagy", 0.85),
    ("KEGG_Mitophagy", "mitophagy", 0.85),
    ("KEGG_mTOR", "upstream_regulators", 0.70),
    ("Reactome_mTOR", "upstream_regulators", 0.70),
]

# Tier priority for when multiple match: pick highest modifier
_TIER_PRIORITY = {
    "core_machinery": 1.0,
    "lysosomal": 1.0,
    "mitophagy": 0.85,
    "docking": 0.85,
    "upstream_regulators": 0.70,
    "default": 0.80,
}


def _map_db_sources_to_tier(db_sources: str) -> str:
    """Infer best scoring tier from pipe-delimited Database Sources string."""
    if not db_sources or db_sources == "nan":
        return "default"
    matched = set()
    for keyword, tier, _ in _SOURCE_TIER_RULES:
        if keyword in db_sources:
            matched.add(tier)
    if not matched:
        return "default"
    return max(matched, key=lambda t: _TIER_PRIORITY.get(t, 0))


# ------------------------------------------------------------------ #
# Load supplementary gene list from Excel
# ------------------------------------------------------------------ #

def load_supplementary_genes(raw_dir: Path, list_cfg: dict) -> dict:
    """Load gene list, build category map, return structured dict."""
    filepath = raw_dir / list_cfg["file"]
    gene_col = list_cfg["gene_symbol_column"]
    hadb_cat_col = list_cfg["hadb_category_column"]

    # -- All Extra Candidates sheet --
    all_df = pd.read_excel(
        filepath, sheet_name=list_cfg["sheets"]["all_candidates"]
    )

    # Expand slash-separated gene symbols
    all_genes = []
    gene_to_db_sources = {}
    for _, row in all_df.iterrows():
        raw = str(row[gene_col]).strip()
        db_src = str(row.get("Database Sources", ""))
        if "/" in raw:
            for g in raw.split("/"):
                g = g.strip()
                all_genes.append(g)
                gene_to_db_sources[g] = db_src
        else:
            all_genes.append(raw)
            gene_to_db_sources[raw] = db_src

    # -- HADb Core sheet --
    hadb_df = pd.read_excel(
        filepath, sheet_name=list_cfg["sheets"]["hadb_core"]
    )
    hadb_categories = {}
    for _, row in hadb_df.iterrows():
        gene = str(row[gene_col]).strip()
        cat = str(row[hadb_cat_col]).strip()
        hadb_categories[gene] = cat

    # -- Build category map: HADb first, then DB sources, fallback default --
    category_map = {}
    for gene in all_genes:
        if gene in hadb_categories:
            category_map[gene] = _map_hadb_category(hadb_categories[gene])
        elif gene in gene_to_db_sources:
            category_map[gene] = _map_db_sources_to_tier(
                gene_to_db_sources[gene]
            )
        else:
            category_map[gene] = "default"

    # De-duplicate while preserving order
    seen = set()
    unique_genes = []
    for g in all_genes:
        if g not in seen:
            unique_genes.append(g)
            seen.add(g)

    return {
        "genes": unique_genes,
        "category_map": category_map,
        "hadb_categories": hadb_categories,
        "label": list_cfg.get("label", "supplementary"),
        "description": list_cfg.get("description", ""),
    }


# ------------------------------------------------------------------ #
# Rescore supplementary genes
# ------------------------------------------------------------------ #

def rescore_genes(
    evidence: pd.DataFrame,
    genes: list,
    category_map: dict,
    category_modifiers: dict,
    weights: dict,
) -> tuple[pd.DataFrame, list]:
    """Recalculate scores for supplementary genes.

    Returns (rescored_df, not_detected_list).
    """
    ev_genes = set(evidence["human_gene_symbol"])
    found = [g for g in genes if g in ev_genes]
    not_detected = [g for g in genes if g not in ev_genes]

    if not found:
        return pd.DataFrame(), not_detected

    mask = evidence["human_gene_symbol"].isin(found)
    df = evidence[mask].copy()

    # Skip genes already in R1 — they already have autophagy scores
    already_r1 = df[df["in_r1"] == True]["human_gene_symbol"].tolist()
    if already_r1:
        logger.info(f"  Skipping {len(already_r1)} genes already in R1: {already_r1}")
        df = df[df["in_r1"] != True].copy()

    if df.empty:
        return pd.DataFrame(), not_detected

    # Store originals
    df["original_score_autophagy"] = df["score_autophagy"]
    df["original_composite_score"] = df["composite_score"]

    # Recalculate score_autophagy
    def new_autophagy(gene):
        tier = category_map.get(gene, "default")
        mod = category_modifiers.get(tier, category_modifiers.get("default", 0.80))
        return 1.0 * mod

    df["score_autophagy"] = df["human_gene_symbol"].apply(new_autophagy)
    df["r1b_tier"] = df["human_gene_symbol"].map(category_map)

    # Recompute composite score
    def recompute(row):
        components = {
            "mouse_csf": row["score_mouse_csf"],
            "human_csf": row["score_human_csf"],
            "ev": row["score_ev"],
            "brain_plausibility": row["score_brain"],
            "autophagy_membership": row["score_autophagy"],
        }
        return compute_composite_score(components, weights, row["penalties"])

    df["composite_score"] = df.apply(recompute, axis=1)
    df["score_delta"] = df["composite_score"] - df["original_composite_score"]
    df["in_r1b"] = True

    df = df.sort_values("composite_score", ascending=False).reset_index(drop=True)
    return df, not_detected


# ------------------------------------------------------------------ #
# Apply hard gates
# ------------------------------------------------------------------ #

def apply_gates(df: pd.DataFrame, config: dict) -> pd.DataFrame:
    """Apply 4 core panel hard gates to R1b candidates."""
    panel_cfg = config["core_panel"]
    mask = pd.Series(True, index=df.index)

    req_human = panel_cfg.get("require_human_csf_tier", ["A", "B"])
    mask &= df["d11_tier"].isin(req_human)
    n_g1 = mask.sum()

    req_mouse = panel_cfg.get("require_mouse_csf_tier", ["A", "B"])
    mask &= df["mouse_csf_best_tier"].isin(req_mouse)
    n_g2 = mask.sum()

    # Gate 3: in_r1b (always True for rescored genes, but explicit)
    mask &= df["in_r1b"] == True
    n_g3 = mask.sum()

    mask &= df["likely_plasma_derived"] == False
    n_g4 = mask.sum()

    logger.info(f"  Gate cascade (R1b):")
    logger.info(f"    Input:               {len(df)}")
    logger.info(f"    G1 (D11 tier A/B):   {n_g1}")
    logger.info(f"    G2 (mouse tier A/B): {n_g2}")
    logger.info(f"    G3 (in_r1b):         {n_g3}")
    logger.info(f"    G4 (not plasma):     {n_g4}")

    return df[mask].sort_values("composite_score", ascending=False)


# ------------------------------------------------------------------ #
# Combined ranking
# ------------------------------------------------------------------ #

def build_combined_ranking(
    r1_core: pd.DataFrame,
    r1b_qualified: pd.DataFrame,
) -> pd.DataFrame:
    """Merge R1 core candidates + R1b qualified into a single ranking."""
    orig = r1_core.copy()
    orig["source"] = "R1"
    if "in_r1b" not in orig.columns:
        orig["in_r1b"] = False

    new = r1b_qualified.copy()
    new["source"] = "R1b"

    # Use only columns present in both
    common_cols = list(set(orig.columns) & set(new.columns))
    combined = pd.concat(
        [orig[common_cols], new[common_cols]], ignore_index=True
    )
    combined = combined.sort_values(
        "composite_score", ascending=False
    ).reset_index(drop=True)
    combined["combined_rank"] = range(1, len(combined) + 1)

    return combined


# ------------------------------------------------------------------ #
# Report writer
# ------------------------------------------------------------------ #

def write_report(
    rescored: pd.DataFrame,
    qualified: pd.DataFrame,
    not_detected: list,
    combined: pd.DataFrame,
    category_map: dict,
    label: str,
    description: str,
    output_path: Path,
):
    """Write human-readable summary report."""
    n_total = len(rescored) + len(not_detected)
    lines = [
        "=" * 70,
        f"SUPPLEMENTARY GENE SCORING REPORT",
        f"List: {label}",
        f"  {description}",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "=" * 70,
        "",
        "-" * 50,
        "INPUT SUMMARY",
        "-" * 50,
        f"  Total supplementary genes:     {n_total}",
        f"  Found in pipeline (rescored):  {len(rescored)}",
        f"  Not detected in any dataset:   {len(not_detected)}",
        "",
    ]

    # Tier distribution
    if len(rescored) > 0:
        tier_counts = rescored["r1b_tier"].value_counts()
        lines += ["-" * 50, "ASSIGNED TIER DISTRIBUTION (found genes)", "-" * 50]
        for tier, count in tier_counts.items():
            mod = _TIER_PRIORITY.get(tier, 0.80)
            lines.append(f"  {tier:<25s} {count:>3d}  (modifier {mod:.2f})")
        lines.append("")

    # Not detected
    if not_detected:
        lines += ["-" * 50, "GENES NOT DETECTED IN ANY DATASET", "-" * 50]
        for g in sorted(not_detected):
            tier = category_map.get(g, "default")
            lines.append(f"  {g:<15s} (tier: {tier})")
        lines.append("")

    # Rescored table
    if len(rescored) > 0:
        lines += [
            "-" * 50,
            "RESCORED RESULTS (sorted by new composite score)",
            "-" * 50,
            f"{'Gene':<14s} {'New':>7s} {'Old':>7s} {'Delta':>7s} "
            f"{'Tier':<20s} {'D11':>4s} {'Mouse':>6s} {'EV':>4s} {'Brain':>6s}",
            "-" * 75,
        ]
        for _, row in rescored.iterrows():
            lines.append(
                f"{row['human_gene_symbol']:<14s} "
                f"{row['composite_score']:>7.3f} "
                f"{row['original_composite_score']:>7.3f} "
                f"{row['score_delta']:>+7.3f} "
                f"{row['r1b_tier']:<20s} "
                f"{row['d11_tier']:>4s} "
                f"{row['mouse_csf_best_tier']:>6s} "
                f"{'Y' if row['ev_present'] else 'N':>4s} "
                f"{'Y' if row['d6_brain_detected'] else 'N':>6s}"
            )
        lines.append("")

    # Core qualified
    lines += [
        "-" * 50,
        f"CORE PANEL QUALIFIED ({len(qualified)} R1b genes pass all 4 gates)",
        "-" * 50,
    ]
    if len(qualified) > 0:
        for _, row in qualified.iterrows():
            lines.append(
                f"  {row['human_gene_symbol']:<14s} "
                f"score={row['composite_score']:.3f}  "
                f"D11={row['d11_tier']}  "
                f"mouse={row['mouse_csf_best_tier']}  "
                f"tier={row['r1b_tier']}"
            )
    else:
        lines.append("  (none)")
    lines.append("")

    # Combined ranking excerpt
    if len(combined) > 0:
        n_r1 = (combined["source"] == "R1").sum()
        n_r1b = (combined["source"] == "R1b").sum()
        lines += [
            "-" * 50,
            f"COMBINED RANKING: {n_r1} R1 + {n_r1b} R1b = {len(combined)} total",
            "-" * 50,
            "",
            "Top 30:",
        ]
        for _, row in combined.head(30).iterrows():
            marker = " *R1b*" if row["source"] == "R1b" else ""
            lines.append(
                f"  {row['combined_rank']:>3d}. {row['human_gene_symbol']:<14s} "
                f"score={row['composite_score']:.3f}  [{row['source']}]{marker}"
            )
        lines.append("")

        # Where do R1b genes land?
        r1b_rows = combined[combined["source"] == "R1b"]
        if len(r1b_rows) > 0:
            lines += ["R1b genes in combined ranking:"]
            for _, row in r1b_rows.iterrows():
                lines.append(
                    f"  rank {row['combined_rank']:>3d}: "
                    f"{row['human_gene_symbol']:<14s} "
                    f"score={row['composite_score']:.3f}"
                )
            lines.append("")

    lines += ["=" * 70, "END OF REPORT", "=" * 70, ""]
    output_path.write_text("\n".join(lines))


# ------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------ #

def main():
    logger.info("=" * 60)
    logger.info("STEP 10: Supplementary Gene Scoring (R1b)")
    logger.info("=" * 60)

    config = load_config(SCRIPTS_DIR / "config.yaml")
    raw_dir = resolve_path(SCRIPTS_DIR, config["paths"]["raw_dir"])
    derived_dir = resolve_path(SCRIPTS_DIR, config["paths"]["derived_dir"])
    output_dir = resolve_path(SCRIPTS_DIR, config["paths"]["output_dir"])
    log_dir = resolve_path(SCRIPTS_DIR, config["paths"]["log_dir"])

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    fh = logging.FileHandler(log_dir / f"10_supplementary_scoring_{timestamp}.log")
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logging.getLogger().addHandler(fh)

    scoring_cfg = config["scoring"]
    weights = scoring_cfg["weights"]
    category_modifiers = scoring_cfg["autophagy_category_modifiers"]

    # Load existing evidence scores (all 9,561 proteins)
    evidence_path = derived_dir / "evidence_scores" / "evidence_scores.tsv"
    evidence = pd.read_csv(evidence_path, sep="\t")
    logger.info(f"Loaded {len(evidence)} proteins from evidence_scores")

    # Load existing R1 core candidates (for combined ranking)
    candidates_path = output_dir / "candidates_ranked.tsv"
    candidates = pd.read_csv(candidates_path, sep="\t")

    # Build the current 122 R1 core: pass all 4 gates
    panel_cfg = config["core_panel"]
    r1_mask = (
        candidates["d11_tier"].isin(panel_cfg["require_human_csf_tier"])
        & candidates["mouse_csf_best_tier"].isin(panel_cfg["require_mouse_csf_tier"])
        & (candidates["in_r1"] == True)
        & (candidates["likely_plasma_derived"] == False)
    )
    r1_core = candidates[r1_mask].sort_values(
        "composite_score", ascending=False
    ).reset_index(drop=True)
    logger.info(f"Loaded {len(r1_core)} R1 core candidates")

    # Process each supplementary list
    supp_lists = config.get("supplementary_gene_lists", {})
    if not supp_lists:
        logger.warning("No supplementary_gene_lists in config — nothing to do")
        return

    for list_id, list_cfg in supp_lists.items():
        logger.info(f"\n{'=' * 50}")
        logger.info(f"Processing: {list_id}")
        logger.info(f"  Label: {list_cfg.get('label', '')}")
        logger.info(f"{'=' * 50}")

        # 1. Load gene list
        supp = load_supplementary_genes(raw_dir, list_cfg)
        logger.info(f"  {len(supp['genes'])} unique genes loaded")
        logger.info(f"  {len(supp['hadb_categories'])} HADb Core genes")

        tier_dist = pd.Series(supp["category_map"]).value_counts()
        for tier, n in tier_dist.items():
            logger.info(f"    {tier}: {n}")

        # 2. Rescore
        rescored, not_detected = rescore_genes(
            evidence, supp["genes"], supp["category_map"],
            category_modifiers, weights,
        )
        logger.info(f"  Rescored: {len(rescored)} | Not detected: {len(not_detected)}")

        # 3. Apply hard gates
        qualified = pd.DataFrame()
        if len(rescored) > 0:
            qualified = apply_gates(rescored, config)
            logger.info(f"  Core qualified: {len(qualified)}")

        # 4. Combined ranking
        combined = build_combined_ranking(r1_core, qualified)
        n_r1b = (combined["source"] == "R1b").sum()
        logger.info(f"  Combined: {len(combined)} ({len(r1_core)} R1 + {n_r1b} R1b)")

        # 5. Sanitise label for filenames
        safe = list_id.lower().replace(" ", "_").replace("/", "_")

        # 6. Save outputs
        if len(rescored) > 0:
            rescored.to_csv(
                output_dir / f"supplementary_{safe}_rescored.tsv",
                sep="\t", index=False,
            )
        if len(qualified) > 0:
            qualified.to_csv(
                output_dir / f"supplementary_{safe}_core_qualified.tsv",
                sep="\t", index=False,
            )
        if not_detected:
            nd_df = pd.DataFrame({
                "gene_symbol": not_detected,
                "assigned_tier": [supp["category_map"].get(g, "default") for g in not_detected],
                "status": "not_detected_in_any_dataset",
            })
            nd_df.to_csv(
                output_dir / f"supplementary_{safe}_not_detected.tsv",
                sep="\t", index=False,
            )
        if len(combined) > 0:
            combined.to_csv(
                output_dir / f"supplementary_{safe}_combined_ranking.tsv",
                sep="\t", index=False,
            )

        # 7. Write report
        write_report(
            rescored, qualified, not_detected, combined,
            supp["category_map"], supp["label"], supp["description"],
            output_dir / f"supplementary_{safe}_report.txt",
        )

        logger.info(f"  Outputs saved to {output_dir}/supplementary_{safe}_*")

    logger.info("\n[10_supplementary_gene_scoring] DONE")


if __name__ == "__main__":
    main()
