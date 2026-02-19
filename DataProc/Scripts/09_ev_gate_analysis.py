#!/usr/bin/env python3
"""Step 9: EV hard-gate sensitivity analysis.

Tests what happens if EV support (ev_present) is made a hard pass/fail
requirement for core panel membership, instead of just a scoring weight.

This is a READ-ONLY analysis — it loads existing pipeline outputs and
re-applies gating logic in-memory. No existing files are modified.

Outputs:
- Outputs/ev_gate_shortlist.tsv         (top 80 under EV hard-gate)
- Outputs/ev_gate_comparison.txt        (diff report vs current shortlist)
"""

import logging
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
logger = logging.getLogger("09_ev_gate_analysis")


def main():
    logger.info("=" * 60)
    logger.info("STEP 9: EV Hard-Gate Sensitivity Analysis")
    logger.info("=" * 60)

    # ── Load config & paths ──
    config = load_config(SCRIPTS_DIR / "config.yaml")
    output_dir = resolve_path(SCRIPTS_DIR, config["paths"]["output_dir"])

    # ── Load existing pipeline outputs ──
    candidates = pd.read_csv(output_dir / "candidates_ranked.tsv", sep="\t")
    current_shortlist = pd.read_csv(
        output_dir / "core_panel_shortlist.tsv", sep="\t"
    )
    logger.info(f"Loaded {len(candidates)} ranked candidates")
    logger.info(f"Loaded {len(current_shortlist)} current shortlist proteins")

    current_genes = set(current_shortlist["human_gene_symbol"])

    # ── Apply all 5 gates (4 current + EV) ──
    panel_cfg = config["core_panel"]
    required_human_tiers = panel_cfg.get("require_human_csf_tier", ["A", "B"])
    required_mouse_tiers = panel_cfg.get("require_mouse_csf_tier", ["A", "B"])
    max_size = panel_cfg.get("max_size", 80)

    mask = pd.Series(True, index=candidates.index)

    # Gate 1: Human CSF tier
    mask &= candidates["d11_tier"].isin(required_human_tiers)
    n_after_g1 = mask.sum()

    # Gate 2: Mouse CSF tier
    mask &= candidates["mouse_csf_best_tier"].isin(required_mouse_tiers)
    n_after_g2 = mask.sum()

    # Gate 3: R1 membership
    mask &= candidates["in_r1"] == True
    n_after_g3 = mask.sum()

    # Gate 4: Plasma exclusion
    mask &= candidates["likely_plasma_derived"] == False
    n_after_g4 = mask.sum()

    # Gate 5: EV requirement (NEW)
    mask &= candidates["ev_present"] == True
    n_after_g5 = mask.sum()

    logger.info("\nGate cascade:")
    logger.info(f"  After G1 (human CSF tier A/B):  {n_after_g1}")
    logger.info(f"  After G2 (mouse CSF tier A/B):  {n_after_g2}")
    logger.info(f"  After G3 (in_r1=True):          {n_after_g3}")
    logger.info(f"  After G4 (not plasma):           {n_after_g4}")
    logger.info(f"  After G5 (ev_present=True):      {n_after_g5}")

    # ── Build EV-gated shortlist ──
    ev_candidates = candidates[mask].sort_values(
        "composite_score", ascending=False
    )
    ev_shortlist = ev_candidates.head(max_size).copy()
    ev_shortlist = ev_shortlist.reset_index(drop=True)

    ev_genes = set(ev_shortlist["human_gene_symbol"])

    # ── Comparison ──
    dropped = current_genes - ev_genes
    gained = ev_genes - current_genes
    retained = current_genes & ev_genes

    logger.info(f"\n{'─' * 50}")
    logger.info("COMPARISON: Current top 80 vs EV-gated top 80")
    logger.info(f"{'─' * 50}")
    logger.info(f"  Retained:  {len(retained)}")
    logger.info(f"  Dropped:   {len(dropped)}")
    logger.info(f"  Gained:    {len(gained)}")

    # ── Build detailed comparison report ──
    lines = [
        "=" * 70,
        "EV HARD-GATE SENSITIVITY ANALYSIS",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "=" * 70,
        "",
        "QUESTION: What if ev_present=True is a hard gate for core panel?",
        "",
        "─" * 50,
        "GATE CASCADE (cumulative filtering)",
        "─" * 50,
        f"  All candidates:                    {len(candidates)}",
        f"  After G1 (human CSF tier A/B):     {n_after_g1}",
        f"  After G2 (mouse CSF tier A/B):     {n_after_g2}",
        f"  After G3 (in_r1=True):             {n_after_g3}",
        f"  After G4 (not plasma-derived):     {n_after_g4}  ← current pipeline",
        f"  After G5 (ev_present=True):        {n_after_g5}  ← with EV gate",
        f"  Shortlist (top {max_size}):               {len(ev_shortlist)}",
        "",
        f"  Proteins lost by adding EV gate:   {n_after_g4 - n_after_g5} "
        f"({n_after_g4} → {n_after_g5})",
        "",
        "─" * 50,
        "SUMMARY",
        "─" * 50,
        f"  Retained in top {max_size}:   {len(retained)} proteins",
        f"  Dropped from top {max_size}:  {len(dropped)} proteins",
        f"  Gained into top {max_size}:   {len(gained)} proteins",
        "",
    ]

    # Dropped proteins detail
    if dropped:
        lines += [
            "─" * 50,
            f"PROTEINS DROPPED FROM TOP {max_size} (lost EV gate)",
            "─" * 50,
            f"{'Gene':<12} {'Curr Rank':>10} {'Score':>7} {'Category':<25} "
            f"{'EV':>3} {'Incl':>4}",
            "-" * 70,
        ]
        dropped_df = current_shortlist[
            current_shortlist["human_gene_symbol"].isin(dropped)
        ].sort_values("composite_score", ascending=False)

        for _, row in dropped_df.iterrows():
            cat = str(row.get("r1_category", ""))[:23]
            ev = "Y" if row.get("ev_present", False) else "N"
            ic = int(row.get("inclusion_count", 0))
            lines.append(
                f"{row['human_gene_symbol']:<12} "
                f"{int(row['rank']):>10} "
                f"{row['composite_score']:>7.3f} "
                f"{cat:<25} "
                f"{ev:>3} "
                f"{ic:>4}"
            )
        lines.append("")

    # Gained proteins detail
    if gained:
        lines += [
            "─" * 50,
            f"PROTEINS GAINED INTO TOP {max_size} (replace dropped)",
            "─" * 50,
            f"{'Gene':<12} {'New Rank':>10} {'Score':>7} {'Category':<25} "
            f"{'EV':>3} {'Incl':>4}",
            "-" * 70,
        ]
        gained_df = ev_shortlist[
            ev_shortlist["human_gene_symbol"].isin(gained)
        ].sort_values("composite_score", ascending=False)

        for _, row in gained_df.iterrows():
            cat = str(row.get("r1_category", ""))[:23]
            ev = "Y" if row.get("ev_present", False) else "N"
            ic = int(row.get("inclusion_count", 0))
            lines.append(
                f"{row['human_gene_symbol']:<12} "
                f"{int(row['rank']):>10} "
                f"{row['composite_score']:>7.3f} "
                f"{cat:<25} "
                f"{ev:>3} "
                f"{ic:>4}"
            )
        lines.append("")

    # Rank shifts for retained proteins
    lines += [
        "─" * 50,
        "RANK SHIFTS (retained proteins, sorted by current rank)",
        "─" * 50,
        f"{'Gene':<12} {'Curr':>5} {'New':>5} {'Shift':>6} {'Score':>7}",
        "-" * 45,
    ]

    # Build rank lookup for both lists
    current_rank = dict(
        zip(
            current_shortlist["human_gene_symbol"],
            range(1, len(current_shortlist) + 1),
        )
    )
    new_rank = dict(
        zip(
            ev_shortlist["human_gene_symbol"],
            range(1, len(ev_shortlist) + 1),
        )
    )

    shifts = []
    for gene in retained:
        cr = current_rank[gene]
        nr = new_rank[gene]
        score = ev_shortlist.loc[
            ev_shortlist["human_gene_symbol"] == gene, "composite_score"
        ].iloc[0]
        shifts.append((gene, cr, nr, nr - cr, score))

    shifts.sort(key=lambda x: x[1])  # sort by current rank

    for gene, cr, nr, shift, score in shifts:
        arrow = "↑" if shift < 0 else ("↓" if shift > 0 else "=")
        shift_str = f"{arrow}{abs(shift)}" if shift != 0 else "  ="
        lines.append(
            f"{gene:<12} {cr:>5} {nr:>5} {shift_str:>6} {score:>7.3f}"
        )

    # AD model cross-check overlay (if available)
    ad_path = output_dir / "ad_model_crosscheck.tsv"
    if ad_path.exists():
        ad_df = pd.read_csv(ad_path, sep="\t")
        lines += [
            "",
            "─" * 50,
            "AD MODEL OVERLAP (dropped proteins)",
            "─" * 50,
        ]
        dropped_in_ad = ad_df[ad_df["human_gene_symbol"].isin(dropped)]
        if len(dropped_in_ad) > 0:
            for _, row in dropped_in_ad.iterrows():
                n_sig = int(row.get("n_comparisons_significant", 0))
                lines.append(
                    f"  {row['human_gene_symbol']:<12} "
                    f"sig_in={n_sig}/3 AD model comparisons"
                )
                if n_sig > 0:
                    lines.append(
                        f"    ⚠  This protein shows significant AD-related "
                        f"changes — dropping it may lose disease relevance"
                    )
        else:
            lines.append("  None of the dropped proteins are in the AD model data")

    lines += ["", "=" * 70, "END OF REPORT", "=" * 70, ""]

    # ── Save outputs ──
    ev_shortlist.to_csv(
        output_dir / "ev_gate_shortlist.tsv", sep="\t", index=False
    )
    logger.info(f"\nSaved EV-gated shortlist: {output_dir / 'ev_gate_shortlist.tsv'}")

    report_text = "\n".join(lines)
    (output_dir / "ev_gate_comparison.txt").write_text(report_text)
    logger.info(f"Saved comparison report: {output_dir / 'ev_gate_comparison.txt'}")

    # Print report to stdout
    print("\n" + report_text)

    logger.info("\n✓ Step 9 complete.")


if __name__ == "__main__":
    main()
