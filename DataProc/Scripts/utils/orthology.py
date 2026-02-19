"""Orthology mapping utilities: gprofiler, mygene, caching."""

import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ------------------------------------------------------------------ #
# gprofiler orthology mapping
# ------------------------------------------------------------------ #

def map_mouse_to_human_gprofiler(gene_symbols: list[str],
                                   organism_from: str = "mmusculus",
                                   organism_to: str = "hsapiens",
                                   batch_size: int = 500,
                                   max_retries: int = 3) -> pd.DataFrame:
    """Map gene symbols between species using gprofiler g:Orth endpoint.

    Parameters
    ----------
    gene_symbols : list of gene symbols to map
    organism_from : source organism (default: mmusculus)
    organism_to : target organism (default: hsapiens)
    batch_size : genes per API call
    max_retries : retry attempts on failure

    Returns
    -------
    DataFrame with columns: input_symbol, output_symbol, orthology_type,
    n_orthologs, description
    """
    from gprofiler import GProfiler

    gp = GProfiler(return_dataframe=True)
    all_results = []

    # Process in batches
    batches = [gene_symbols[i:i + batch_size]
               for i in range(0, len(gene_symbols), batch_size)]

    logger.info(f"Mapping {len(gene_symbols)} genes via gprofiler "
                f"({organism_from} -> {organism_to}), {len(batches)} batch(es)")

    for batch_idx, batch in enumerate(batches):
        for attempt in range(max_retries):
            try:
                result = gp.orth(
                    organism=organism_from,
                    query=batch,
                    target=organism_to,
                )
                all_results.append(result)
                logger.info(f"  Batch {batch_idx + 1}/{len(batches)}: "
                           f"{len(result)} mappings returned")
                break
            except Exception as e:
                wait = 2 ** (attempt + 1)
                logger.warning(f"  gprofiler batch {batch_idx + 1} attempt {attempt + 1} "
                             f"failed: {e}. Retrying in {wait}s...")
                time.sleep(wait)
                if attempt == max_retries - 1:
                    logger.error(f"  gprofiler batch {batch_idx + 1} failed after "
                               f"{max_retries} attempts")
                    raise

    if not all_results:
        return pd.DataFrame(columns=["input_symbol", "output_symbol",
                                     "orthology_type", "n_orthologs", "description"])

    combined = pd.concat(all_results, ignore_index=True)

    # gprofiler returns: incoming (input mouse symbol), converted (mouse Ensembl ID),
    # ortholog_ensg (human Ensembl ID), name (human gene symbol)
    # We want output_symbol = human gene symbol (from 'name' column)
    col_map = {
        "incoming": "input_symbol",
        "name": "output_symbol",
    }
    result_df = combined.rename(columns={k: v for k, v in col_map.items() if k in combined.columns})
    # Drop the mouse Ensembl ID column ('converted') to avoid confusion
    if "converted" in result_df.columns:
        result_df = result_df.drop(columns=["converted"])

    # Determine orthology type per input gene
    if "input_symbol" in result_df.columns:
        # Count orthologs per input
        ortho_counts = result_df.groupby("input_symbol")["output_symbol"].nunique()
        result_df["n_orthologs"] = result_df["input_symbol"].map(ortho_counts)

        # Classify: 1:1, 1:many
        result_df["orthology_type"] = result_df["n_orthologs"].apply(
            lambda n: "1:1" if n == 1 else "1:many" if n > 1 else "none"
        )

    # Handle N/A or "N/A" output symbols (unmapped genes)
    na_mask = (
        result_df["output_symbol"].isna()
        | (result_df["output_symbol"].astype(str).str.upper() == "N/A")
        | (result_df["output_symbol"].astype(str).str.strip() == "")
        | (result_df["output_symbol"].astype(str) == "None")
    )
    result_df.loc[na_mask, "orthology_type"] = "none"
    result_df.loc[na_mask, "n_orthologs"] = 0

    logger.info(f"  Total mappings: {len(result_df)}, "
                f"unique inputs: {result_df['input_symbol'].nunique()}")

    return result_df


# ------------------------------------------------------------------ #
# mygene cross-validation
# ------------------------------------------------------------------ #

def cross_validate_mygene(gene_pairs: list[tuple],
                            source_species: str = "mouse",
                            target_species: str = "human",
                            batch_size: int = 1000) -> pd.DataFrame:
    """Cross-validate orthology mappings using mygene.info.

    Parameters
    ----------
    gene_pairs : list of (source_symbol, candidate_target_symbol) tuples
    source_species : source species name
    target_species : target species name

    Returns
    -------
    DataFrame with: input_symbol, candidate_symbol, mygene_confirmed (bool),
    mygene_alternatives (list)
    """
    import mygene

    mg = mygene.MyGeneInfo()
    source_symbols = list(set(p[0] for p in gene_pairs))
    species_tax = {"mouse": 10090, "human": 9606}

    logger.info(f"Cross-validating {len(source_symbols)} genes via mygene")

    results = []
    for i in range(0, len(source_symbols), batch_size):
        batch = source_symbols[i:i + batch_size]
        try:
            query_results = mg.querymany(
                batch,
                scopes="symbol",
                fields="symbol,homologene",
                species=species_tax.get(source_species, source_species),
                returnall=True,
            )
            for hit in query_results.get("out", []):
                symbol = hit.get("query", "")
                homologene = hit.get("homologene", {})
                homolog_genes = []
                if isinstance(homologene, dict):
                    for gene_entry in homologene.get("genes", []):
                        if isinstance(gene_entry, list) and len(gene_entry) >= 2:
                            tax_id, gene_id = gene_entry[0], gene_entry[1]
                            if tax_id == species_tax.get(target_species, 0):
                                homolog_genes.append(str(gene_id))

                results.append({
                    "input_symbol": symbol,
                    "mygene_homologs": homolog_genes,
                })
        except Exception as e:
            logger.warning(f"  mygene batch failed: {e}")

    # Match against candidate pairs
    mygene_map = {r["input_symbol"]: r["mygene_homologs"] for r in results}
    validation = []
    for source, candidate in gene_pairs:
        homologs = mygene_map.get(source, [])
        validation.append({
            "input_symbol": source,
            "candidate_symbol": candidate,
            "mygene_confirmed": candidate in homologs or candidate.upper() in [h.upper() for h in homologs],
            "mygene_alternatives": ";".join(homologs),
        })

    return pd.DataFrame(validation)


# ------------------------------------------------------------------ #
# Cache management
# ------------------------------------------------------------------ #

def load_orthology_cache(cache_path: Path) -> pd.DataFrame:
    """Load cached orthology mappings from TSV."""
    if cache_path.exists():
        df = pd.read_csv(cache_path, sep="\t")
        logger.info(f"Loaded orthology cache: {len(df)} entries from {cache_path}")
        return df
    return pd.DataFrame(columns=["input_symbol", "output_symbol",
                                 "orthology_type", "n_orthologs"])


def save_orthology_cache(mappings: pd.DataFrame, cache_path: Path) -> None:
    """Save orthology mappings to TSV cache."""
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    mappings.to_csv(cache_path, sep="\t", index=False)
    logger.info(f"Saved orthology cache: {len(mappings)} entries to {cache_path}")


# ------------------------------------------------------------------ #
# Orthology resolution
# ------------------------------------------------------------------ #

def resolve_orthology(gene_symbol: str, species: str,
                       cache: pd.DataFrame,
                       policy: str = "expand") -> list[dict]:
    """Look up orthology mapping and apply resolution policy.

    Parameters
    ----------
    gene_symbol : gene symbol to look up
    species : species of the gene ('mouse' or 'human')
    cache : orthology cache DataFrame
    policy : 'expand' (all orthologs) or 'exclude' (skip ambiguous)

    Returns
    -------
    List of dicts with: human_gene_symbol, orthology_ambiguous, orthology_type
    """
    if species == "human":
        return [{"human_gene_symbol": gene_symbol,
                 "orthology_ambiguous": False,
                 "orthology_type": "same_species"}]

    if cache.empty:
        return [{"human_gene_symbol": None,
                 "orthology_ambiguous": False,
                 "orthology_type": "unmapped"}]

    matches = cache[cache["input_symbol"].str.upper() == gene_symbol.upper()]

    if matches.empty:
        return [{"human_gene_symbol": None,
                 "orthology_ambiguous": False,
                 "orthology_type": "unmapped"}]

    # Filter out N/A outputs
    valid = matches[
        matches["output_symbol"].notna()
        & (matches["output_symbol"].astype(str).str.upper() != "N/A")
        & (matches["output_symbol"].astype(str).str.strip() != "")
    ]

    if valid.empty:
        return [{"human_gene_symbol": None,
                 "orthology_ambiguous": False,
                 "orthology_type": "unmapped"}]

    ortho_type = valid.iloc[0].get("orthology_type", "unknown")

    if ortho_type == "1:1" or len(valid) == 1:
        return [{"human_gene_symbol": valid.iloc[0]["output_symbol"],
                 "orthology_ambiguous": False,
                 "orthology_type": ortho_type}]

    # 1:many case
    if policy == "expand":
        return [
            {"human_gene_symbol": row["output_symbol"],
             "orthology_ambiguous": True,
             "orthology_type": "1:many"}
            for _, row in valid.iterrows()
        ]
    elif policy == "exclude":
        return [{"human_gene_symbol": None,
                 "orthology_ambiguous": True,
                 "orthology_type": "1:many_excluded"}]
    else:
        raise ValueError(f"Unknown orthology policy: {policy}")


def precurate_r1_orthology(r1_genes: set, cache: pd.DataFrame) -> pd.DataFrame:
    """Identify R1 autophagy genes with non-trivial orthology mappings.

    Returns DataFrame of genes with 1:many or many:1 mappings.
    """
    records = []
    for gene in sorted(r1_genes):
        matches = cache[cache["output_symbol"].str.upper() == gene.upper()]
        if matches.empty:
            # Try as input
            matches = cache[cache["input_symbol"].str.upper() == gene.upper()]

        if matches.empty:
            continue

        n_orthologs = int(matches.iloc[0].get("n_orthologs", 1))
        if n_orthologs > 1:
            records.append({
                "human_gene_symbol": gene,
                "n_orthologs": n_orthologs,
                "orthology_type": matches.iloc[0].get("orthology_type", "unknown"),
                "all_orthologs": ";".join(matches["input_symbol"].unique())
                    if "input_symbol" in matches.columns else "",
            })

    result = pd.DataFrame(records)
    logger.info(f"R1 orthology pre-curation: {len(result)} genes with non-trivial mappings")
    return result
