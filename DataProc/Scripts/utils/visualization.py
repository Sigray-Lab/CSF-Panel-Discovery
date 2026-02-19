"""Visualization utilities: UpSet plots, heatmaps, Venn diagrams, etc."""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns

logger = logging.getLogger(__name__)


def plot_upset(sets_dict: dict, title: str, output_path: Path,
               fig_format: str = "pdf") -> None:
    """Generate UpSet plot for multi-set intersections."""
    from upsetplot import from_memberships, UpSet

    # Build membership list
    all_elements = set()
    for s in sets_dict.values():
        all_elements.update(s)

    memberships = []
    for element in all_elements:
        member_of = [name for name, s in sets_dict.items() if element in s]
        if member_of:
            memberships.append(member_of)

    if not memberships:
        logger.warning("No data for UpSet plot")
        return

    data = from_memberships(memberships)
    upset = UpSet(data, show_counts=True, sort_by="cardinality")

    fig = plt.figure(figsize=(12, 6))
    upset.plot(fig=fig)
    plt.suptitle(title, fontsize=14, y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close()
    logger.info(f"  Saved UpSet plot: {output_path}")


def plot_detection_heatmap(detection_matrix: pd.DataFrame,
                             output_path: Path,
                             tier_column_prefix: str = "",
                             max_proteins: int = 100) -> None:
    """Heatmap of panel candidates x datasets.

    detection_matrix: rows = proteins, columns = dataset tiers
    Color by tier: A=3, B=2, C=1, absent=0
    """
    tier_map = {"A": 3, "B": 2, "C": 1, "absent": 0}

    # Limit to top N proteins
    plot_df = detection_matrix.head(max_proteins).copy()

    # Convert tiers to numeric
    for col in plot_df.columns:
        if col != "human_gene_symbol":
            plot_df[col] = plot_df[col].map(tier_map).fillna(0)

    plot_df = plot_df.set_index("human_gene_symbol")

    fig, ax = plt.subplots(figsize=(max(8, len(plot_df.columns) * 0.8),
                                     max(6, len(plot_df) * 0.25)))
    cmap = sns.color_palette(["#f0f0f0", "#fdd49e", "#fdbb84", "#e34a33"], as_cmap=True)
    sns.heatmap(plot_df, cmap=cmap, vmin=0, vmax=3,
                cbar_kws={"ticks": [0, 1, 2, 3],
                          "label": "Tier (0=absent, 1=C, 2=B, 3=A)"},
                ax=ax, linewidths=0.5)
    ax.set_title("Detection Tier Heatmap (Panel Candidates x Datasets)")
    ax.set_ylabel("Protein")
    ax.set_xlabel("Dataset")
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close()
    logger.info(f"  Saved detection heatmap: {output_path}")


def plot_score_distribution(scores_autophagy: pd.Series,
                              scores_non_autophagy: pd.Series,
                              output_path: Path) -> None:
    """Overlapping histograms of evidence scores: autophagy vs non-autophagy."""
    fig, ax = plt.subplots(figsize=(8, 5))

    bins = np.linspace(0, 1, 40)
    ax.hist(scores_non_autophagy, bins=bins, alpha=0.6, color="steelblue",
            label=f"Non-autophagy (n={len(scores_non_autophagy)})", density=True)
    ax.hist(scores_autophagy, bins=bins, alpha=0.6, color="firebrick",
            label=f"Autophagy/R1 (n={len(scores_autophagy)})", density=True)

    ax.set_xlabel("Composite Evidence Score")
    ax.set_ylabel("Density")
    ax.set_title("Evidence Score Distribution")
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close()
    logger.info(f"  Saved score distribution: {output_path}")


def plot_sensitivity_stability(rank_matrix: pd.DataFrame,
                                 output_path: Path,
                                 top_n: int = 50) -> None:
    """Candidate rank consistency across parameter grid.

    rank_matrix: columns = parameter combinations, rows = proteins,
    values = ranks.
    """
    plot_df = rank_matrix.head(top_n)

    fig, ax = plt.subplots(figsize=(10, max(6, top_n * 0.25)))
    cmap = sns.color_palette("YlOrRd_r", as_cmap=True)
    sns.heatmap(plot_df, cmap=cmap, ax=ax, linewidths=0.3,
                cbar_kws={"label": "Rank"})
    ax.set_title(f"Rank Stability Across Parameter Grid (Top {top_n})")
    ax.set_ylabel("Protein")
    ax.set_xlabel("Parameter Combination")
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close()
    logger.info(f"  Saved sensitivity stability: {output_path}")


def plot_module_enrichment(enrichment_df: pd.DataFrame,
                             output_path: Path) -> None:
    """Bar chart of R1 enrichment per module."""
    if enrichment_df.empty:
        logger.warning("No enrichment data for module plot")
        return

    sig = enrichment_df[enrichment_df["pvalue_adj"] < 0.05].copy()
    if sig.empty:
        sig = enrichment_df.head(10).copy()

    sig = sig.sort_values("fold_enrichment", ascending=True)

    fig, ax = plt.subplots(figsize=(8, max(4, len(sig) * 0.4)))
    colors = ["firebrick" if p < 0.05 else "steelblue"
              for p in sig["pvalue_adj"]]
    ax.barh(sig["module_id"].astype(str), sig["fold_enrichment"], color=colors)
    ax.set_xlabel("Fold Enrichment (R1 Autophagy Genes)")
    ax.set_ylabel("Module")
    ax.set_title("Module Enrichment for Autophagy/Lysosome Genes")
    ax.axvline(x=1.0, color="gray", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close()
    logger.info(f"  Saved module enrichment: {output_path}")


def plot_comparison_venn(new_panel: set, r2: set, previous_329: set,
                           output_path: Path) -> None:
    """3-way Venn diagram: new panel vs R2 vs previous 329."""
    from matplotlib_venn import venn3

    fig, ax = plt.subplots(figsize=(8, 6))
    venn3(
        [new_panel, r2, previous_329],
        set_labels=("New Panel", "R2 (Antibody)", "Previous 329"),
        ax=ax,
    )
    ax.set_title("Panel Comparison")
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close()
    logger.info(f"  Saved comparison Venn: {output_path}")


def plot_null_distribution(observed: int, null_counts: list,
                             output_path: Path) -> None:
    """Histogram of null simulation convergence counts with observed marked."""
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.hist(null_counts, bins=30, alpha=0.7, color="steelblue",
            label="Null simulations")
    ax.axvline(x=observed, color="firebrick", linewidth=2, linestyle="--",
               label=f"Observed ({observed})")

    # p-value annotation
    n_greater = sum(1 for c in null_counts if c >= observed)
    p_value = n_greater / max(len(null_counts), 1)
    ax.text(0.95, 0.95, f"p = {p_value:.4f}\nFold = {observed / max(np.mean(null_counts), 1):.1f}x",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=11, bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    ax.set_xlabel("Number of Convergent Genes")
    ax.set_ylabel("Count (Simulations)")
    ax.set_title("Chance Expectation: Autophagy Gene Convergence")
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close()
    logger.info(f"  Saved null distribution: {output_path}")
