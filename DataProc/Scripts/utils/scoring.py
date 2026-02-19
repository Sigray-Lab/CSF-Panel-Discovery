"""Evidence score computation and R1 category mapping."""

import logging

import numpy as np

logger = logging.getLogger(__name__)


# ------------------------------------------------------------------ #
# R1 category mapping
# ------------------------------------------------------------------ #

# Keyword-based mapping from heterogeneous R1 category strings
# to the 5 scoring tiers defined in config.yaml
CATEGORY_KEYWORDS = {
    "core_machinery": [
        "ULK1", "ULK complex", "ATG12", "ATG8", "ATG9", "ATG16",
        "ATG5", "ATG7", "ATG3", "ATG10", "ATG14",
        "Beclin", "PI3K complex", "VPS34", "VPS15",
        "Cargo receptor", "cargo receptor",
        "Ubiquitin-like conjugation",
        "Phagophore", "phagophore",
        "Autophagosome formation", "autophagosome formation",
        "Nucleation", "nucleation",
        "Elongation", "elongation",
        "Closure", "closure",
    ],
    "lysosomal": [
        "Lysosom", "lysosom",
        "acid hydrolase", "Acid hydrolase",
        "cathepsin", "Cathepsin",
        "ATP hydrolysis", "V-type", "v-type", "V-ATPase", "v-ATPase",
        "lysosomal acidification", "Lysosomal acidification",
        "Lysosome Biogenesis", "lysosome biogenesis",
        "Lysosome membrane", "lysosome membrane",
        "CLN", "Neuronal ceroid",
        "Sphingolipid", "sphingolipid",
        "Glycosidase", "glycosidase",
        "Lipase", "lipase",
        "Protease", "protease",
        "Sulfatase", "sulfatase",
        "Phospholipase", "phospholipase",
        "GCase", "Glucocerebrosidase",
    ],
    "mitophagy": [
        "Mitophagy", "mitophagy",
        "PINK1", "Parkin", "BNIP3", "NIX", "FUNDC1",
        "mitochondrial quality",
    ],
    "docking": [
        "SNARE", "Tether", "tether",
        "Autophagosome-lysosome fusion", "autophagosome-lysosome fusion",
        "HOPS", "Rab protein", "RAB",
        "Adaptor", "adaptor",
        "Docking", "docking",
        "Fusion", "fusion",
    ],
    "upstream_regulators": [
        "mTOR", "MTOR", "AMPK",
        "Positive regulator", "Negative regulator",
        "positive regulator", "negative regulator",
        "Transcription Factor", "transcription factor",
        "TFEB", "TFE3",
        "Candidate kinase", "candidate kinase",
        "PGC1", "Signaling", "signaling",
        "upstream", "Upstream",
        "PI3K/AKT", "Insulin",
        "Nutrient sensing", "nutrient sensing",
        "autophagic amino acid",
    ],
}


def categorize_r1_category(raw_category: str) -> str:
    """Map a heterogeneous R1 category string to a scoring tier.

    Returns one of: core_machinery, lysosomal, mitophagy, docking,
    upstream_regulators, default
    """
    if not raw_category or str(raw_category).strip() in ("", "nan", "None"):
        return "default"

    cat_str = str(raw_category)

    for tier, keywords in CATEGORY_KEYWORDS.items():
        for kw in keywords:
            if kw in cat_str:
                return tier

    return "default"


# ------------------------------------------------------------------ #
# Score component computation
# ------------------------------------------------------------------ #

def compute_mouse_csf_evidence(protein_data: dict,
                                 mouse_csf_datasets: list[str]) -> float:
    """Compute mouse CSF evidence component.

    Formula: (datasets_detected / total_mouse_CSF_datasets)
             * mean(pct_detected across detected datasets)

    D1 contributes binary only (detected=1 or 0).
    Range: 0-1
    """
    total_datasets = len(mouse_csf_datasets)
    if total_datasets == 0:
        return 0.0

    datasets_detected = 0
    pct_values = []

    for ds_id in mouse_csf_datasets:
        ds_info = protein_data.get(ds_id)
        if ds_info is None:
            continue

        detected = ds_info.get("detected", False)
        if detected:
            datasets_detected += 1
            if ds_id == "D1":
                # D1 is presence-only, contribute binary
                pct_values.append(1.0)
            else:
                pct = ds_info.get("pct_detected", 0.0)
                pct_values.append(pct)

    if datasets_detected == 0:
        return 0.0

    breadth = datasets_detected / total_datasets
    mean_pct = sum(pct_values) / len(pct_values) if pct_values else 0.0

    return breadth * mean_pct


def compute_human_csf_evidence(protein_data: dict,
                                 tier_scores: dict,
                                 d10_bonus: float = 0.1) -> float:
    """Compute human CSF evidence component.

    D11 detectability tier -> score (A=1.0, B=0.7, C=0.3, absent=0).
    D10 bonus: +0.1 if detected.
    D12: validation flag only, not scored.
    Range: 0-1.1
    """
    d11_info = protein_data.get("D11", {})
    d11_tier = d11_info.get("detectability_tier", "absent")
    score = tier_scores.get(d11_tier, tier_scores.get("absent", 0.0))

    # D10 bonus
    d10_info = protein_data.get("D10", {})
    if d10_info.get("detected", False):
        score += d10_bonus

    return score


def compute_ev_evidence(gene_symbol: str, ev_data: dict) -> float:
    """Compute EV evidence component.

    ev_present * (1 - rank_percentile)
    Lower rank = higher abundance = higher score.
    Range: 0-1
    """
    if gene_symbol is None:
        return 0.0

    rank_percentiles = ev_data.get("rank_percentiles", {})
    genes = ev_data.get("genes", set())

    gene_upper = str(gene_symbol).upper()
    # Check presence (case-insensitive)
    if gene_upper not in {g.upper() for g in genes}:
        return 0.0

    # Find rank percentile (case-insensitive lookup)
    for g, pct in rank_percentiles.items():
        if g.upper() == gene_upper:
            return 1.0 - pct

    # Present but no rank info
    return 0.5  # default mid-range if present but rank unknown


def compute_brain_plausibility(protein_data: dict) -> float:
    """D6 brain lysate: detected = 1.0, absent = 0."""
    d6_info = protein_data.get("D6", {})
    return 1.0 if d6_info.get("detected", False) else 0.0


def compute_autophagy_membership(gene_symbol: str,
                                   r1_genes: set,
                                   category_map: dict,
                                   category_modifiers: dict) -> float:
    """Compute autophagy membership component.

    In R1 = base 1.0 * category modifier.
    Not in R1 = 0.
    Range: 0-1
    """
    if gene_symbol is None:
        return 0.0

    gene_upper = str(gene_symbol).upper()
    r1_upper = {g.upper() for g in r1_genes}

    if gene_upper not in r1_upper:
        return 0.0

    # Find original-case gene for category lookup
    raw_category = None
    for g in r1_genes:
        if g.upper() == gene_upper:
            raw_category = category_map.get(g)
            break

    if raw_category is None:
        raw_category = category_map.get(gene_symbol, "")

    tier = categorize_r1_category(raw_category)
    modifier = category_modifiers.get(tier, category_modifiers.get("default", 0.80))

    return 1.0 * modifier


def compute_penalties(is_ambiguous: bool, orthology_ambiguous: bool,
                       likely_plasma: bool, penalty_config: dict) -> float:
    """Sum applicable penalties."""
    total = 0.0
    if is_ambiguous:
        total += penalty_config.get("ambiguous_group_unresolved", -0.10)
    if orthology_ambiguous:
        total += penalty_config.get("orthology_ambiguous_unresolved", -0.10)
    if likely_plasma:
        total += penalty_config.get("likely_plasma_derived", -0.20)
    return total


def compute_composite_score(components: dict, weights: dict,
                              penalties: float) -> float:
    """Weighted sum of components + penalties, clamped to [0, 1]."""
    score = 0.0
    for key, value in components.items():
        w = weights.get(key, 0.0)
        score += w * value

    score += penalties
    return max(0.0, min(1.0, score))
