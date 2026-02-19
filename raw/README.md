# Raw Data Sources

This folder contains the source proteomics data required to run the pipeline. These files are **not included** in the repository (see licensing/consent notes below). To reproduce the analysis, obtain the files listed below and place them in this directory with the exact folder structure shown.

## Required files

### Reference lists

| File | Description | Source |
|------|-------------|--------|
| `Full list with autophagy proteins curated.xlsx` | Curated autophagy/lysosome gene list (R1, 599 genes) | Study authors |
| `WT vs KO EV protein list from highest to lowest abundance.xlsx` | EV secretome reference (SH-SY5Y cells) | Study authors |

### Pooled datasets

| File | Sheets used | Description | Source |
|------|-------------|-------------|--------|
| `Pooled datasets blasting overview.xlsx` | `mouse_csf_PXD053698` (D2), `mouse_csf_PXD053568_5xFAD ` (D3, note trailing space), `mouse_isf_PDX048864` (D9), `human_csf_PXD053698` (D10) | Multi-study pooled proteomics | ProteomeXchange: PXD053698, PXD053568, PXD048864 |

### ProteomeXchange datasets

| Folder / File | Dataset ID | Species | Biofluid | Role |
|---------------|-----------|---------|----------|------|
| `The list of all identified proteins from lable-free MS of mouse CSF(Per).xlsx` | D1 | Mouse | CSF | Mouse CSF (label-free) |
| `Project PXD032782/APPPS1_BI_2_2Pep.xlsx` | D4 (PXD032782) | Mouse | CSF | Mouse CSF (APP/PS1) |
| `Project PXD065438/Output_CSF-mk_Trem2/CSF-mk_Trem2.pg_matrix.tsv` | D5 (PXD065438) | Mouse | CSF | Mouse CSF (TREM2) |
| `Project PXD065438/Output_Brain-mk_Trem2/AGCH2501_brain.pg_matrix.tsv` | D6 (PXD065438) | Mouse | Brain | Brain plausibility |
| `Project PXD052590/report.pg_matrix.tsv` | D7 (PXD052590) | Mouse | CSF | Mouse CSF |
| `Project PXD048864/report.pg_matrix_CSF.tsv` | D8 (PXD048864) | Mouse | CSF | Mouse CSF |
| `Project Astral/Astral protome_discovery and repliction_DACboth013_49D_1_7_3_Sc01.tsv` | D11 | Human | CSF | Human CSF primary |
| `Project Astral/Astral proteome_repliction only_DACRP014_49D_1_7_2_Sc03.tsv` | D12 | Human | CSF | Human CSF validation |

### AD model cross-check (Step 08)

| File | Description |
|------|-------------|
| `Table S3 Alterations in proteins identifed in all samples of AppNL-G-F NL-G-F vs Appwt wt.xlsx` | App knock-in AD mouse model CSF proteomics (3 comparisons) |

### Supplementary gene lists (Step 10)

| File | Description |
|------|-------------|
| `Misc/novel_candidates_INCLUDE_56_extra_from_Edison_Cosmos_run.xlsx` | 56 novel + 9 HADb core autophagy genes |
| `Misc/CSF_Panel_Pathway_Additions.xlsx` | Pathway-specific candidate proteins (sEH, cGAS-STING, Sigma-1R, tau degradation) |
| `Misc/TDP43_additional_biomarkers.xlsx` | TDP-43-related biomarkers (STMN2, GRN) |

## Data access

Most mouse CSF datasets are publicly available via [ProteomeXchange](https://www.proteomexchange.org/) using the accession numbers listed above. The human Astral CSF data and curated reference lists should be requested from the study authors.

## Notes

- File names must match exactly (including spaces and capitalisation)
- The D3 sheet name has a **trailing space**: `"mouse_csf_PXD053568_5xFAD "`
- The D9 sheet uses `PDX` instead of `PXD`: `"mouse_isf_PDX048864"`
- See `DataProc/Scripts/config.yaml` for the complete dataset registry with all column mappings
