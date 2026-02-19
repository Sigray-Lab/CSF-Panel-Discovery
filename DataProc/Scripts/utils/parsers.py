"""Dataset-specific parsers for all proteomics data formats.

Handles 5 distinct formats:
- maxquant_fasta (D1): FASTA header GN= regex parsing
- maxquant_preparsed (D2/D3/D10): pre-parsed gene_symbol columns
- diann_tsv (D5/D6/D7/D8): DIA-NN pg_matrix TSV files
- diann_excel (D4/D9): DIA-NN format in Excel workbooks
- diann_astral (D11/D12): Astral DIA-NN with extra columns
- ev_rank (EV): rank-ordered gene symbol lists
"""

import re
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

logger = logging.getLogger(__name__)


# ------------------------------------------------------------------ #
# Config loading
# ------------------------------------------------------------------ #

def load_config(config_path: str | Path) -> dict:
    """Load and return config.yaml as dict."""
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    with open(config_path) as f:
        config = yaml.safe_load(f)
    return config


def resolve_path(raw_dir: str | Path, relative_path: str) -> Path:
    """Resolve a dataset file path relative to raw_dir."""
    return Path(raw_dir) / relative_path


# ------------------------------------------------------------------ #
# Sample column identification
# ------------------------------------------------------------------ #

def identify_sample_columns(df: pd.DataFrame, ds_config: dict) -> list[str]:
    """Identify intensity/sample columns by excluding known meta columns.

    Strategy varies by format:
    - maxquant_fasta/maxquant_preparsed: columns starting with intensity_prefix
    - diann_tsv/diann_astral: columns not in meta_columns list
    - diann_excel: columns starting with intensity_prefix, skipping extra prefixes
    """
    fmt = ds_config["format"]
    all_cols = list(df.columns)

    # Columns added by the pipeline (parsers, QC, orthology, standardisation)
    # that should never be treated as sample/intensity columns
    _PIPELINE_COLUMNS = {
        "group_id", "uniprot_ids", "group_members", "gene_symbol",
        "is_ambiguous_group", "human_gene_symbol", "orthology_ambiguous",
        "orthology_type", "n_detected", "pct_detected", "median_intensity",
        "above_intensity_floor", "detectability_tier", "likely_plasma_derived",
        "dataset_id", "species", "protein",
    }

    if fmt in ("maxquant_fasta", "maxquant_preparsed"):
        prefix = ds_config["intensity_prefix"]
        return [c for c in all_cols if c.startswith(prefix)]

    if fmt == "diann_excel":
        meta = set(ds_config.get("meta_columns", []))
        exclude = meta | _PIPELINE_COLUMNS
        skip_prefixes = ds_config.get("extra_column_prefixes_to_skip", [])
        prefix = ds_config.get("intensity_prefix")
        sample_cols = []
        for c in all_cols:
            if c in exclude:
                continue
            if any(c.startswith(sp) for sp in skip_prefixes):
                continue
            if prefix and not c.startswith(prefix):
                continue
            sample_cols.append(c)
        # If no prefix filter, take everything not in meta/skip/pipeline
        if not prefix:
            sample_cols = [c for c in all_cols
                           if c not in exclude
                           and not any(c.startswith(sp) for sp in skip_prefixes)]
        return sample_cols

    if fmt in ("diann_tsv", "diann_astral"):
        meta = set(ds_config.get("meta_columns", []))
        exclude = meta | _PIPELINE_COLUMNS
        return [c for c in all_cols if c not in exclude]

    raise ValueError(f"Unknown format for sample column detection: {fmt}")


# ------------------------------------------------------------------ #
# Protein group ambiguity detection
# ------------------------------------------------------------------ #

def _parse_gene_list(genes_str: str | float) -> list[str]:
    """Split a semicolon-delimited gene string into unique gene symbols."""
    if pd.isna(genes_str) or str(genes_str).strip() == "":
        return []
    parts = [g.strip() for g in str(genes_str).split(";") if g.strip()]
    return parts


def _check_ambiguous_group(genes_str: str | float) -> bool:
    """True if the gene list contains >1 unique gene symbol."""
    genes = _parse_gene_list(genes_str)
    unique = set(g.upper() for g in genes)
    return len(unique) > 1


def _first_gene(genes_str: str | float) -> str | None:
    """Extract first gene symbol from semicolon-delimited list."""
    genes = _parse_gene_list(genes_str)
    return genes[0] if genes else None


# ------------------------------------------------------------------ #
# MaxQuant parsers
# ------------------------------------------------------------------ #

def parse_maxquant_fasta(filepath: Path, sheet: str, ds_config: dict) -> pd.DataFrame:
    """Parse D1-format MaxQuant files requiring FASTA header GN= regex.

    Returns DataFrame with standardised columns + per-sample intensities.
    """
    logger.info(f"Parsing MaxQuant FASTA format: {filepath} sheet={sheet}")
    df = pd.read_excel(filepath, sheet_name=sheet, engine="openpyxl")
    logger.info(f"  Raw shape: {df.shape}")

    # Extract gene symbols from Fasta headers via GN= regex
    def extract_gene(fasta_str):
        if pd.isna(fasta_str):
            return None
        # Take first FASTA entry if semicolon-separated
        first_entry = str(fasta_str).split(";")[0]
        match = re.search(r"GN=(\w+)", first_entry)
        return match.group(1) if match else None

    df["gene_symbol"] = df[ds_config["id_columns"]["fasta_headers"]].apply(extract_gene)

    # Parse protein IDs: extract first UniProt accession from sp|ACC|NAME format
    def extract_accessions(pid_str):
        if pd.isna(pid_str):
            return []
        accessions = []
        for entry in str(pid_str).split(";"):
            entry = entry.strip()
            # sp|P12345|NAME or tr|P12345|NAME format
            match = re.match(r"(?:sp|tr)\|([A-Z0-9_-]+)\|", entry)
            if match:
                accessions.append(match.group(1))
            elif entry and not entry.startswith("REV__") and not entry.startswith("CON__"):
                accessions.append(entry)
        return accessions

    pid_col = ds_config["id_columns"]["protein_ids"]
    df["_accessions"] = df[pid_col].apply(extract_accessions)
    df["group_id"] = df[pid_col].astype(str)
    df["uniprot_ids"] = df["_accessions"].apply(lambda x: ";".join(x))
    df["group_members"] = df["uniprot_ids"]

    # Gene-level ambiguity: extract all genes from all FASTA headers
    def extract_all_genes(fasta_str):
        if pd.isna(fasta_str):
            return []
        return re.findall(r"GN=(\w+)", str(fasta_str))

    fasta_col = ds_config["id_columns"]["fasta_headers"]
    df["_all_genes"] = df[fasta_col].apply(extract_all_genes)
    df["is_ambiguous_group"] = df["_all_genes"].apply(
        lambda genes: len(set(g.upper() for g in genes)) > 1
    )

    # Identify sample columns
    sample_cols = identify_sample_columns(df, ds_config)
    logger.info(f"  Found {len(sample_cols)} sample columns")

    # Replace zeros with NaN in intensity columns
    for col in sample_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
        df[col] = df[col].replace(0, np.nan)

    # Clean up temp columns
    df = df.drop(columns=["_accessions", "_all_genes"], errors="ignore")

    return df


def parse_maxquant_preparsed(filepath: Path, sheet: str, ds_config: dict) -> pd.DataFrame:
    """Parse D2/D3/D10-format pre-parsed MaxQuant sheets.

    These already have gene_symbol columns. The uniprot_id column often
    contains truncated/invalid values and should not be relied upon.
    """
    logger.info(f"Parsing MaxQuant preparsed: {filepath} sheet='{sheet}'")
    df = pd.read_excel(filepath, sheet_name=sheet, engine="openpyxl")
    logger.info(f"  Raw shape: {df.shape}")

    gene_col = ds_config["id_columns"]["gene_symbol"]
    df["gene_symbol"] = df[gene_col].astype(str).str.strip()
    df.loc[df["gene_symbol"].isin(["nan", "", "None"]), "gene_symbol"] = np.nan

    # Use gene_symbol as group_id since uniprot_id may be unreliable
    df["group_id"] = df["gene_symbol"]
    df["uniprot_ids"] = ""
    df["group_members"] = df["gene_symbol"]
    df["is_ambiguous_group"] = False  # one gene per row in preparsed format

    # Identify sample columns
    sample_cols = identify_sample_columns(df, ds_config)
    logger.info(f"  Found {len(sample_cols)} sample columns")

    # Replace zeros with NaN
    for col in sample_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
        df[col] = df[col].replace(0, np.nan)

    return df


# ------------------------------------------------------------------ #
# DIA-NN parsers
# ------------------------------------------------------------------ #

def _parse_diann_ids(df: pd.DataFrame, ds_config: dict) -> pd.DataFrame:
    """Parse DIA-NN-style identifiers: Protein.Group, Genes, optional Protein.Ids."""
    id_cols = ds_config["id_columns"]

    # Protein.Group -> group_id
    pg_col = id_cols.get("protein_group", "Protein.Group")
    df["group_id"] = df[pg_col].astype(str)

    # Protein.Ids (optional) -> uniprot_ids
    if "protein_ids" in id_cols and id_cols["protein_ids"] in df.columns:
        df["uniprot_ids"] = df[id_cols["protein_ids"]].astype(str)
    else:
        df["uniprot_ids"] = df["group_id"]

    df["group_members"] = df["group_id"]

    # Genes -> gene_symbol (first gene) + ambiguity check
    genes_col = id_cols.get("genes", "Genes")
    if genes_col in df.columns:
        df["gene_symbol"] = df[genes_col].apply(_first_gene)
        df["is_ambiguous_group"] = df[genes_col].apply(_check_ambiguous_group)
    else:
        df["gene_symbol"] = np.nan
        df["is_ambiguous_group"] = False
        logger.warning(f"  No Genes column found ({genes_col})")

    return df


def parse_diann_tsv(filepath: Path, ds_config: dict) -> pd.DataFrame:
    """Parse DIA-NN pg_matrix TSV files (D5, D6, D7, D8)."""
    logger.info(f"Parsing DIA-NN TSV: {filepath}")
    df = pd.read_csv(filepath, sep="\t")
    logger.info(f"  Raw shape: {df.shape}")

    df = _parse_diann_ids(df, ds_config)

    sample_cols = identify_sample_columns(df, ds_config)
    logger.info(f"  Found {len(sample_cols)} sample columns")

    # Ensure numeric
    for col in sample_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def parse_diann_excel(filepath: Path, ds_config: dict) -> pd.DataFrame:
    """Parse DIA-NN format from Excel workbooks (D4, D9).

    D4: has LFQ, Precursor, and Peptides column groups - use only LFQ.
    D9: from pooled workbook, has 2 completely empty sample columns.
    """
    sheet = ds_config.get("sheet")
    logger.info(f"Parsing DIA-NN Excel: {filepath} sheet={sheet}")
    df = pd.read_excel(filepath, sheet_name=sheet, engine="openpyxl")
    logger.info(f"  Raw shape: {df.shape}")

    df = _parse_diann_ids(df, ds_config)

    sample_cols = identify_sample_columns(df, ds_config)
    logger.info(f"  Found {len(sample_cols)} sample columns")

    # Ensure numeric
    for col in sample_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def parse_diann_astral(filepath: Path, ds_config: dict) -> pd.DataFrame:
    """Parse Astral DIA-NN format (D11, D12).

    Key differences:
    - Has extra 'protein' column at position 0
    - D12 has 4 precomputed stats columns
    - Log2 intensities (do not replace zeros)
    - Use float32 for memory efficiency
    """
    logger.info(f"Parsing DIA-NN Astral: {filepath}")

    # Read with optimized dtypes
    df = pd.read_csv(filepath, sep="\t", low_memory=False)
    logger.info(f"  Raw shape: {df.shape}")

    df = _parse_diann_ids(df, ds_config)

    sample_cols = identify_sample_columns(df, ds_config)
    logger.info(f"  Found {len(sample_cols)} sample columns")

    # Convert sample columns to float32 for memory efficiency
    for col in sample_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce").astype(np.float32)

    return df


# ------------------------------------------------------------------ #
# EV dataset parser
# ------------------------------------------------------------------ #

def parse_ev_dataset(filepath: Path, ds_config: dict) -> dict:
    """Parse EV dataset with multiple sheets.

    Returns dict with:
    - 'raw_wt': list of gene symbols from WT column (ranked)
    - 'raw_ko': list of gene symbols from KO column (ranked)
    - 'rank_percentiles': dict mapping gene_symbol -> rank_percentile (best of WT/KO)
    - 'precomputed': DataFrame from '330 after blast with full' sheet
    """
    logger.info(f"Parsing EV dataset: {filepath}")
    sheets = ds_config["sheets"]

    # Raw list
    raw_df = pd.read_excel(filepath, sheet_name=sheets["raw"], engine="openpyxl")
    raw_wt = raw_df["WT"].dropna().astype(str).str.strip().tolist()
    raw_ko = raw_df["KO"].dropna().astype(str).str.strip().tolist()
    logger.info(f"  Raw WT genes: {len(raw_wt)}, KO genes: {len(raw_ko)}")

    # Compute rank percentiles from raw lists
    rank_percentiles = {}
    for i, gene in enumerate(raw_wt):
        if gene and gene != "nan":
            rank_percentiles[gene] = i / max(len(raw_wt) - 1, 1)
    for i, gene in enumerate(raw_ko):
        if gene and gene != "nan":
            ko_pct = i / max(len(raw_ko) - 1, 1)
            # Keep best (lowest) rank percentile
            if gene not in rank_percentiles or ko_pct < rank_percentiles[gene]:
                rank_percentiles[gene] = ko_pct

    # Precomputed 329-protein intersection
    precomputed = pd.read_excel(filepath, sheet_name=sheets["precomputed"],
                                engine="openpyxl")
    logger.info(f"  Precomputed intersection: {len(precomputed)} proteins")

    # All EV genes (union of WT and KO)
    all_genes = set(g for g in raw_wt + raw_ko if g and g != "nan")

    return {
        "raw_wt": raw_wt,
        "raw_ko": raw_ko,
        "rank_percentiles": rank_percentiles,
        "genes": all_genes,
        "precomputed": precomputed,
    }


# ------------------------------------------------------------------ #
# Reference input parsers
# ------------------------------------------------------------------ #

def parse_r1_reference(filepath: Path, ref_config: dict) -> dict:
    """Parse R1 autophagy reference.

    Returns dict with:
    - 'genes': set of unique gene symbols from deduplicated sheet
    - 'full_list': DataFrame from full sheet with categories
    - 'category_map': dict mapping gene_symbol -> primary category string
    """
    logger.info(f"Parsing R1 reference: {filepath}")
    cols = ref_config["columns"]

    # Deduplicated sheet: unique gene symbols
    dedup_df = pd.read_excel(filepath, sheet_name=ref_config["deduplicated_sheet"],
                             engine="openpyxl")
    gene_col = cols["gene_symbol"]
    genes = set(dedup_df[gene_col].dropna().astype(str).str.strip())
    genes.discard("")
    genes.discard("nan")
    logger.info(f"  Unique genes from deduplicated sheet: {len(genes)}")

    # Full list with categories
    full_df = pd.read_excel(filepath, sheet_name=ref_config["full_sheet"],
                            engine="openpyxl")
    logger.info(f"  Full list rows: {len(full_df)}")

    # Build category map: gene -> first encountered category
    category_map = {}
    for _, row in full_df.iterrows():
        gene = str(row.get(gene_col, "")).strip()
        cat = str(row.get(cols["category"], "")).strip()
        if gene and gene != "nan" and gene not in category_map:
            category_map[gene] = cat

    return {
        "genes": genes,
        "full_list": full_df,
        "category_map": category_map,
    }


def parse_r2_reference(filepath: Path, ref_config: dict) -> pd.DataFrame:
    """Parse R2 overview reference (antibody targets)."""
    logger.info(f"Parsing R2 reference: {filepath} sheet={ref_config['sheet']}")
    df = pd.read_excel(filepath, sheet_name=ref_config["sheet"], engine="openpyxl")
    cols = ref_config["columns"]

    # Standardise column names
    result = pd.DataFrame()
    result["antibody"] = df[cols["antibody"]].astype(str).str.strip()
    result["gene_symbol"] = df[cols["gene_symbol"]].astype(str).str.strip()
    if cols["abbreviation"] in df.columns:
        result["abbreviation"] = df[cols["abbreviation"]].astype(str).str.strip()

    # Drop rows without gene symbols
    result = result[result["gene_symbol"].notna() & (result["gene_symbol"] != "nan")]
    logger.info(f"  R2 targets: {len(result)}")
    return result


# ------------------------------------------------------------------ #
# Standardisation
# ------------------------------------------------------------------ #

def standardize_columns(df: pd.DataFrame, dataset_id: str, species: str,
                        sample_cols: list[str]) -> pd.DataFrame:
    """Ensure DataFrame has all output contract columns.

    Required columns: group_id, uniprot_ids, gene_symbol, group_members,
    is_ambiguous_group, species, dataset_id.

    Adds placeholder columns for downstream steps if not present:
    human_gene_symbol, orthology_ambiguous, detectability_tier,
    likely_plasma_derived, n_detected, pct_detected, median_intensity,
    above_intensity_floor.
    """
    required = {
        "group_id": "",
        "uniprot_ids": "",
        "gene_symbol": np.nan,
        "group_members": "",
        "is_ambiguous_group": False,
    }

    for col, default in required.items():
        if col not in df.columns:
            df[col] = default

    df["dataset_id"] = dataset_id
    df["species"] = species

    # Downstream placeholders
    for col in ["human_gene_symbol", "orthology_ambiguous", "detectability_tier",
                "likely_plasma_derived", "n_detected", "pct_detected",
                "median_intensity", "above_intensity_floor"]:
        if col not in df.columns:
            df[col] = np.nan

    # Ensure standard column order: contract cols first, then sample cols
    contract_cols = [
        "dataset_id", "species", "group_id", "uniprot_ids", "gene_symbol",
        "group_members", "is_ambiguous_group", "human_gene_symbol",
        "orthology_ambiguous", "detectability_tier", "likely_plasma_derived",
        "n_detected", "pct_detected", "median_intensity", "above_intensity_floor",
    ]
    other_cols = [c for c in df.columns if c not in contract_cols and c not in sample_cols]
    ordered = contract_cols + other_cols + sample_cols
    # Only include columns that exist
    ordered = [c for c in ordered if c in df.columns]
    df = df[ordered]

    return df
