#!/usr/bin/env python3
"""Step 5: Peptide feasibility and cross-species conservation.

For top N candidates by score:
1. Fetch UniProt sequences (human + mouse ortholog)
2. In silico tryptic digest
3. Identify proteotypic peptides
4. Cross-species peptide conservation
5. Assign feasibility tiers (I/II/III)

Outputs:
- DerivedData/peptide_feasibility/peptide_feasibility.tsv
"""

import logging
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS_DIR))

from utils.parsers import load_config
from utils.orthology import load_orthology_cache

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("05_peptide_feasibility")


# ------------------------------------------------------------------ #
# UniProt sequence fetching
# ------------------------------------------------------------------ #

def fetch_uniprot_sequence(gene_symbol: str, species: str,
                            cache_dir: Path,
                            max_retries: int = 3) -> str | None:
    """Fetch protein sequence from UniProt REST API with caching.

    Parameters
    ----------
    gene_symbol : gene symbol to look up
    species : 'human' or 'mouse'
    cache_dir : directory for FASTA cache files

    Returns
    -------
    Protein sequence string, or None if not found
    """
    import requests

    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"{species}_{gene_symbol}.fasta"

    # Check cache first
    if cache_file.exists():
        content = cache_file.read_text().strip()
        if content:
            # Parse FASTA: skip header line
            lines = content.split("\n")
            seq = "".join(l.strip() for l in lines if not l.startswith(">"))
            if seq:
                return seq

    # Fetch from UniProt
    organism = "human" if species == "human" else "mouse"
    tax_id = 9606 if species == "human" else 10090

    url = (
        f"https://rest.uniprot.org/uniprotkb/search?"
        f"query=gene_exact:{gene_symbol}+AND+organism_id:{tax_id}"
        f"+AND+reviewed:true"
        f"&format=fasta&size=1"
    )

    for attempt in range(max_retries):
        try:
            resp = requests.get(url, timeout=30)
            if resp.status_code == 200 and resp.text.strip():
                cache_file.write_text(resp.text)
                lines = resp.text.strip().split("\n")
                seq = "".join(l.strip() for l in lines if not l.startswith(">"))
                return seq if seq else None
            elif resp.status_code == 200:
                # Empty response — gene not found
                return None
            else:
                logger.warning(f"  UniProt {gene_symbol} ({species}): "
                             f"HTTP {resp.status_code}")
        except Exception as e:
            logger.warning(f"  UniProt {gene_symbol} ({species}) attempt "
                         f"{attempt + 1}: {e}")
            time.sleep(2 ** attempt)

    return None


# ------------------------------------------------------------------ #
# Tryptic digest
# ------------------------------------------------------------------ #

def tryptic_digest(sequence: str, missed_cleavages: list[int],
                    min_length: int = 7, max_length: int = 25) -> list[str]:
    """In silico tryptic digest using pyteomics.

    Returns list of unique peptide sequences within length bounds.
    """
    from pyteomics import parser

    all_peptides = set()
    for mc in missed_cleavages:
        peptides = parser.cleave(sequence, "trypsin", missed_cleavages=mc)
        for pep in peptides:
            if min_length <= len(pep) <= max_length:
                all_peptides.add(pep)

    return sorted(all_peptides)


# ------------------------------------------------------------------ #
# Proteotypicity check
# ------------------------------------------------------------------ #

def check_proteotypicity(peptides: list[str], all_sequences: dict) -> list[str]:
    """Filter for proteotypic peptides (unique to one protein).

    Simple approach: check if peptide substring occurs in only one
    protein sequence in the reference set.

    Parameters
    ----------
    peptides : candidate peptide sequences
    all_sequences : dict mapping gene_symbol -> protein sequence
        (from the same species)

    Returns
    -------
    List of proteotypic peptides
    """
    proteotypic = []
    for pep in peptides:
        n_hits = sum(1 for seq in all_sequences.values()
                     if pep in seq)
        if n_hits == 1:
            proteotypic.append(pep)
    return proteotypic


# ------------------------------------------------------------------ #
# Conservation check
# ------------------------------------------------------------------ #

def find_conserved_peptides(human_peptides: list[str],
                              mouse_peptides: list[str],
                              allow_mismatch: int = 0) -> list[str]:
    """Find peptides shared between human and mouse.

    Parameters
    ----------
    human_peptides : human proteotypic peptides
    mouse_peptides : mouse proteotypic peptides
    allow_mismatch : number of allowed mismatches (0 = exact match)

    Returns
    -------
    List of conserved peptide sequences (from human set)
    """
    if allow_mismatch == 0:
        mouse_set = set(mouse_peptides)
        return [p for p in human_peptides if p in mouse_set]
    else:
        # Near-identical: allow N substitutions
        conserved = []
        for hp in human_peptides:
            for mp in mouse_peptides:
                if len(hp) == len(mp):
                    mismatches = sum(1 for a, b in zip(hp, mp) if a != b)
                    if mismatches <= allow_mismatch:
                        conserved.append(hp)
                        break
        return conserved


# ------------------------------------------------------------------ #
# Feasibility tier assignment
# ------------------------------------------------------------------ #

def assign_feasibility_tier(n_conserved: int, n_human: int,
                              tier_config: dict) -> str:
    """Assign peptide feasibility tier (I/II/III)."""
    tier1_min = tier_config.get("I", {}).get("min_conserved_proteotypic", 2)
    tier2_min = tier_config.get("II", {}).get("min_human_proteotypic", 1)

    if n_conserved >= tier1_min:
        return "I"
    elif n_human >= tier2_min:
        return "II"
    else:
        return "III"


# ------------------------------------------------------------------ #
# Main
# ------------------------------------------------------------------ #

def main():
    config_path = SCRIPTS_DIR / "config.yaml"
    config = load_config(config_path)
    derived_dir = (SCRIPTS_DIR / config["paths"]["derived_dir"]).resolve()
    log_dir = (SCRIPTS_DIR / config["paths"]["log_dir"]).resolve()
    cache_dir = derived_dir / "uniprot_cache"
    output_dir = derived_dir / "peptide_feasibility"
    output_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    file_handler = logging.FileHandler(log_dir / f"05_peptide_feasibility_{timestamp}.log")
    file_handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logging.getLogger().addHandler(file_handler)

    pf_config = config["peptide_feasibility"]
    missed_cleavages = pf_config["missed_cleavages"]
    min_len = pf_config["peptide_length_min"]
    max_len = pf_config["peptide_length_max"]
    top_n = pf_config["top_n_candidates"]
    tier_config = pf_config["tiers"]

    logger.info("=" * 70)
    logger.info("Step 5: Peptide Feasibility and Cross-Species Conservation")
    logger.info(f"Top N candidates: {top_n}")
    logger.info("=" * 70)

    # ------------------------------------------------------------------
    # 1. Load core panel candidates
    # ------------------------------------------------------------------
    candidates_path = derived_dir / "autophagy_panel" / "core_panel_candidates.tsv"
    if not candidates_path.exists():
        # Fallback to evidence scores
        candidates_path = derived_dir / "evidence_scores" / "evidence_scores.tsv"

    candidates = pd.read_csv(candidates_path, sep="\t")
    candidates = candidates.head(top_n)
    logger.info(f"Processing {len(candidates)} candidates")

    # ------------------------------------------------------------------
    # 2. Load orthology cache for mouse gene lookups
    # ------------------------------------------------------------------
    ortho_cache = load_orthology_cache(derived_dir / "orthology_cache.tsv")

    # Build human->mouse reverse lookup
    human_to_mouse = {}
    if not ortho_cache.empty and "input_symbol" in ortho_cache.columns:
        for _, row in ortho_cache.iterrows():
            out = str(row.get("output_symbol", "")).strip()
            inp = str(row.get("input_symbol", "")).strip()
            if out and out != "nan" and out != "N/A" and inp and inp != "nan":
                if out not in human_to_mouse:
                    human_to_mouse[out] = inp

    # ------------------------------------------------------------------
    # 3. Process each candidate
    # ------------------------------------------------------------------
    results = []
    # Collect all sequences for proteotypicity check
    human_sequences = {}
    mouse_sequences = {}

    # First pass: fetch all sequences
    logger.info("Fetching UniProt sequences...")
    genes_to_process = candidates["human_gene_symbol"].dropna().unique()

    for gene in genes_to_process:
        # Human sequence
        human_seq = fetch_uniprot_sequence(gene, "human", cache_dir)
        if human_seq:
            human_sequences[gene] = human_seq

        # Mouse ortholog sequence
        mouse_gene = human_to_mouse.get(gene)
        if mouse_gene:
            mouse_seq = fetch_uniprot_sequence(mouse_gene, "mouse", cache_dir)
            if mouse_seq:
                mouse_sequences[mouse_gene] = mouse_seq

        # Rate limiting
        time.sleep(0.2)

    logger.info(f"Fetched {len(human_sequences)} human, {len(mouse_sequences)} mouse sequences")

    # Second pass: digest and assess
    logger.info("Running tryptic digest and conservation analysis...")
    for _, row in candidates.iterrows():
        gene = row.get("human_gene_symbol")
        if pd.isna(gene) or gene not in human_sequences:
            results.append({
                "human_gene_symbol": gene,
                "n_tryptic_peptides_human": 0,
                "n_proteotypic_human": 0,
                "n_tryptic_peptides_mouse": 0,
                "n_proteotypic_mouse": 0,
                "n_conserved_proteotypic": 0,
                "feasibility_tier": "III",
                "best_conserved_peptides": "",
            })
            continue

        # Human digest
        human_seq = human_sequences[gene]
        human_peptides = tryptic_digest(human_seq, missed_cleavages, min_len, max_len)

        # Proteotypicity check (simplified: just check against all fetched sequences)
        human_proteotypic = check_proteotypicity(human_peptides, human_sequences)

        # Mouse digest
        mouse_gene = human_to_mouse.get(gene)
        mouse_peptides = []
        mouse_proteotypic = []
        if mouse_gene and mouse_gene in mouse_sequences:
            mouse_seq = mouse_sequences[mouse_gene]
            mouse_peptides = tryptic_digest(mouse_seq, missed_cleavages, min_len, max_len)
            mouse_proteotypic = check_proteotypicity(mouse_peptides, mouse_sequences)

        # Conservation
        conserved = find_conserved_peptides(human_proteotypic, mouse_proteotypic)

        # Feasibility tier
        tier = assign_feasibility_tier(len(conserved), len(human_proteotypic), tier_config)

        results.append({
            "human_gene_symbol": gene,
            "n_tryptic_peptides_human": len(human_peptides),
            "n_proteotypic_human": len(human_proteotypic),
            "n_tryptic_peptides_mouse": len(mouse_peptides),
            "n_proteotypic_mouse": len(mouse_proteotypic),
            "n_conserved_proteotypic": len(conserved),
            "feasibility_tier": tier,
            "best_conserved_peptides": ";".join(conserved[:5]),
            "mouse_ortholog": mouse_gene or "",
        })

    results_df = pd.DataFrame(results)

    # ------------------------------------------------------------------
    # 4. Save output
    # ------------------------------------------------------------------
    output_path = output_dir / "peptide_feasibility.tsv"
    results_df.to_csv(output_path, sep="\t", index=False)
    logger.info(f"\nPeptide feasibility: {output_path}")

    # Summary
    tier_counts = results_df["feasibility_tier"].value_counts()
    logger.info(f"Feasibility tier distribution: {tier_counts.to_dict()}")
    logger.info(f"  Tier I (conserved proteotypic): {tier_counts.get('I', 0)}")
    logger.info(f"  Tier II (human-only): {tier_counts.get('II', 0)}")
    logger.info(f"  Tier III (no good peptides): {tier_counts.get('III', 0)}")

    logger.info("[05_peptide_feasibility] DONE")


if __name__ == "__main__":
    main()
