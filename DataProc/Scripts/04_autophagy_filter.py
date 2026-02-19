#!/usr/bin/env python3
"""Step 4: Autophagy filtering and biological annotation.

1. Validate R1 against GO:0006914 and Reactome R-HSA-9612973
2. Annotate candidates with biological metadata
3. Run enrichment on convergent-but-not-autophagy proteins
4. Define core panel candidates

Outputs:
- DerivedData/autophagy_panel/annotated_candidates.tsv
- DerivedData/autophagy_panel/core_panel_candidates.tsv
- QC/R1_validation_vs_GO_Reactome.tsv
- QC/convergent_non_autophagy_enrichment.tsv
"""

import logging
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS_DIR))

from utils.parsers import load_config, resolve_path, parse_r1_reference
from utils.scoring import categorize_r1_category

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("04_autophagy_filter")


# ------------------------------------------------------------------ #
# R1 validation
# ------------------------------------------------------------------ #

def validate_r1_against_go_reactome(r1_genes: set, qc_dir: Path) -> pd.DataFrame:
    """Cross-reference R1 against GO autophagy and Reactome autophagy terms.

    Uses gseapy to fetch gene sets from Enrichr libraries.
    """
    logger.info("Validating R1 against GO/Reactome...")

    go_autophagy_genes = set()
    reactome_autophagy_genes = set()

    try:
        import gseapy as gp

        # Fetch GO Biological Process autophagy terms
        try:
            go_lib = gp.get_library("GO_Biological_Process_2023")
            for term, genes in go_lib.items():
                term_lower = term.lower()
                if "autophagy" in term_lower or "autophag" in term_lower:
                    go_autophagy_genes.update(g.upper() for g in genes)
            logger.info(f"  GO autophagy terms: {len(go_autophagy_genes)} genes")
        except Exception as e:
            logger.warning(f"  Could not fetch GO library: {e}")

        # Fetch Reactome autophagy terms
        try:
            reactome_lib = gp.get_library("Reactome_2022")
            for term, genes in reactome_lib.items():
                term_lower = term.lower()
                if "autophagy" in term_lower or "autophag" in term_lower:
                    reactome_autophagy_genes.update(g.upper() for g in genes)
            logger.info(f"  Reactome autophagy terms: {len(reactome_autophagy_genes)} genes")
        except Exception as e:
            logger.warning(f"  Could not fetch Reactome library: {e}")

    except ImportError:
        logger.warning("  gseapy not available for GO/Reactome validation")

    # Build comparison table
    r1_upper = {g.upper() for g in r1_genes}
    all_genes = r1_upper | go_autophagy_genes | reactome_autophagy_genes

    records = []
    for gene in sorted(all_genes):
        in_r1 = gene in r1_upper
        in_go = gene in go_autophagy_genes
        in_reactome = gene in reactome_autophagy_genes

        if in_r1 and (in_go or in_reactome):
            status = "confirmed"
        elif in_r1:
            status = "r1_only"
        elif in_go or in_reactome:
            status = "go_reactome_only"
        else:
            status = "unknown"

        records.append({
            "gene_symbol": gene,
            "in_r1": in_r1,
            "in_go_autophagy": in_go,
            "in_reactome_autophagy": in_reactome,
            "status": status,
        })

    result = pd.DataFrame(records)

    # Summary
    for status in ["confirmed", "r1_only", "go_reactome_only"]:
        n = len(result[result["status"] == status])
        logger.info(f"  {status}: {n} genes")

    # Save
    output_path = qc_dir / "R1_validation_vs_GO_Reactome.tsv"
    result.to_csv(output_path, sep="\t", index=False)
    logger.info(f"  Saved: {output_path}")

    return result


# ------------------------------------------------------------------ #
# Biological annotation
# ------------------------------------------------------------------ #

def annotate_candidates(scores: pd.DataFrame, r1: dict) -> pd.DataFrame:
    """Add biological annotation columns.

    Adds:
    - r1_functional_category: from R1 Full list
    - r1_scoring_tier: mapped tier
    - secretory_autophagy_cargo: known secretory autophagy cargoes
    """
    df = scores.copy()

    # R1 functional category
    df["r1_functional_category"] = df["human_gene_symbol"].map(r1["category_map"])

    # Scoring tier
    df["r1_scoring_tier"] = df["r1_functional_category"].apply(categorize_r1_category)

    # Known secretory autophagy cargoes
    SECRETORY_AUTOPHAGY = {
        "SEC22B", "TRIM16", "IL1B", "IL18", "IL33",
        "HMGB1", "ANXA1", "ANXA2", "FTH1", "FTL",
        "GSDMD",
    }
    df["secretory_autophagy_cargo"] = df["human_gene_symbol"].apply(
        lambda g: str(g).upper() in SECRETORY_AUTOPHAGY if pd.notna(g) else False
    )

    return df


# ------------------------------------------------------------------ #
# Enrichment analysis
# ------------------------------------------------------------------ #

def run_enrichment_analysis(gene_list: list[str], qc_dir: Path) -> pd.DataFrame:
    """Run GO/Reactome enrichment on a gene list via gseapy."""
    logger.info(f"Running enrichment analysis on {len(gene_list)} genes...")

    if len(gene_list) < 5:
        logger.warning("  Too few genes for enrichment analysis")
        return pd.DataFrame()

    try:
        import gseapy as gp

        results = gp.enrich(
            gene_list=gene_list,
            gene_sets=["GO_Biological_Process_2023", "Reactome_2022"],
            organism="human",
            outdir=None,
            no_plot=True,
        )

        if results is not None and hasattr(results, "results") and not results.results.empty:
            enrichment = results.results[
                ["Term", "Gene_set", "P-value", "Adjusted P-value",
                 "Overlap", "Genes"]
            ].copy()
            enrichment = enrichment.sort_values("Adjusted P-value")
            logger.info(f"  Found {len(enrichment)} enriched terms "
                       f"({(enrichment['Adjusted P-value'] < 0.05).sum()} significant)")
            return enrichment
        else:
            logger.info("  No significant enrichment found")
            return pd.DataFrame()

    except ImportError:
        logger.warning("  gseapy not available for enrichment analysis")
        return pd.DataFrame()
    except Exception as e:
        logger.warning(f"  Enrichment analysis failed: {e}")
        return pd.DataFrame()


# ------------------------------------------------------------------ #
# Core panel definition
# ------------------------------------------------------------------ #

def define_core_panel(scores: pd.DataFrame, config: dict) -> pd.DataFrame:
    """Apply core panel selection criteria.

    Criteria from config:
    - Human CSF tier A or B in D11
    - Mouse CSF tier A or B in >= 1 dataset
    - In R1 autophagy list
    - Not likely_plasma_derived
    - EV is never required
    """
    panel_cfg = config["core_panel"]

    mask = pd.Series(True, index=scores.index)

    # Human CSF tier requirement
    required_tiers = panel_cfg.get("require_human_csf_tier", ["A", "B"])
    mask &= scores["d11_tier"].isin(required_tiers)

    # Mouse CSF tier requirement
    required_mouse_tiers = panel_cfg.get("require_mouse_csf_tier", ["A", "B"])
    mask &= scores["mouse_csf_best_tier"].isin(required_mouse_tiers)

    # R1 membership
    if panel_cfg.get("require_r1", True):
        mask &= scores["in_r1"] == True

    # Plasma exclusion
    if panel_cfg.get("exclude_plasma", True):
        mask &= scores["likely_plasma_derived"] == False

    core = scores[mask].sort_values("composite_score", ascending=False)
    logger.info(f"Core panel: {len(core)} candidates (criteria: "
                f"D11 tier {required_tiers}, mouse tier {required_mouse_tiers}, "
                f"R1={'required' if panel_cfg.get('require_r1') else 'optional'}, "
                f"plasma={'excluded' if panel_cfg.get('exclude_plasma') else 'included'})")

    return core


def main():
    config_path = SCRIPTS_DIR / "config.yaml"
    config = load_config(config_path)
    raw_dir = (SCRIPTS_DIR / config["paths"]["raw_dir"]).resolve()
    derived_dir = (SCRIPTS_DIR / config["paths"]["derived_dir"]).resolve()
    qc_dir = (SCRIPTS_DIR / config["paths"]["qc_dir"]).resolve()
    output_dir = (SCRIPTS_DIR / config["paths"]["output_dir"]).resolve()
    log_dir = (SCRIPTS_DIR / config["paths"]["log_dir"]).resolve()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    file_handler = logging.FileHandler(log_dir / f"04_autophagy_filter_{timestamp}.log")
    file_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logging.getLogger().addHandler(file_handler)

    panel_dir = derived_dir / "autophagy_panel"
    panel_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 70)
    logger.info("Step 4: Autophagy Filtering and Biological Annotation")
    logger.info("=" * 70)

    # ------------------------------------------------------------------
    # 1. Load scored data and R1
    # ------------------------------------------------------------------
    scores_path = derived_dir / "evidence_scores" / "evidence_scores.tsv"
    scores = pd.read_csv(scores_path, sep="\t")
    logger.info(f"Loaded {len(scores)} scored proteins")

    r1_ref = config["references"]["R1"]
    r1_path = resolve_path(raw_dir, r1_ref["file"])
    r1 = parse_r1_reference(r1_path, r1_ref)

    # ------------------------------------------------------------------
    # 2. R1 validation against GO/Reactome
    # ------------------------------------------------------------------
    validate_r1_against_go_reactome(r1["genes"], qc_dir)

    # ------------------------------------------------------------------
    # 3. Annotate candidates
    # ------------------------------------------------------------------
    annotated = annotate_candidates(scores, r1)
    logger.info(f"Annotated {len(annotated)} proteins")

    # ------------------------------------------------------------------
    # 4. Enrichment of convergent-but-not-autophagy proteins
    # ------------------------------------------------------------------
    convergent_non_r1 = annotated[
        (annotated["n_mouse_csf_datasets"] > 0)
        & (annotated["d11_tier"].isin(["A", "B"]))
        & (~annotated["in_r1"])
    ]
    logger.info(f"Convergent non-R1 proteins: {len(convergent_non_r1)}")

    if len(convergent_non_r1) >= 5:
        enrichment = run_enrichment_analysis(
            convergent_non_r1["human_gene_symbol"].dropna().tolist(),
            qc_dir,
        )
        if not enrichment.empty:
            enr_path = qc_dir / "convergent_non_autophagy_enrichment.tsv"
            enrichment.to_csv(enr_path, sep="\t", index=False)
            logger.info(f"Enrichment results: {enr_path}")
    else:
        logger.info("Too few convergent non-R1 proteins for enrichment")

    # ------------------------------------------------------------------
    # 5. Define core panel
    # ------------------------------------------------------------------
    core_panel = define_core_panel(annotated, config)

    # ------------------------------------------------------------------
    # 6. Save outputs
    # ------------------------------------------------------------------
    annotated_path = panel_dir / "annotated_candidates.tsv"
    annotated.to_csv(annotated_path, sep="\t", index=False)
    logger.info(f"Annotated candidates: {annotated_path}")

    core_path = panel_dir / "core_panel_candidates.tsv"
    core_panel.to_csv(core_path, sep="\t", index=False)
    logger.info(f"Core panel candidates: {core_path} ({len(core_panel)} proteins)")

    # Summary
    logger.info("\nTop 20 core panel candidates:")
    for _, row in core_panel.head(20).iterrows():
        logger.info(f"  {row['human_gene_symbol']:15s} "
                    f"score={row['composite_score']:.4f} "
                    f"D11={row['d11_tier']} "
                    f"mouse={row['n_mouse_csf_datasets']}ds "
                    f"cat={row.get('r1_scoring_tier', 'N/A')}")

    logger.info("[04_autophagy_filter] DONE")


if __name__ == "__main__":
    main()
