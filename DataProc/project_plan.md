# CSF Autophagy/Lysosome Panel Discovery — Pipeline Specification v3 (FINAL)

This document is the **self-sufficient implementation spec** for a reproducible pipeline that identifies and ranks **CSF-detectable markers of autophagy/lysosomal function** with **rodent↔human translation** in mind.

**Key design corrections vs the original plan:**
- **EV evidence is NOT a hard gate** — supportive annotation only.
- Replace Boolean "detected/not detected" with **dataset-aware detectability tiers (A/B/C)** + **weighted evidence scoring**.
- **Accession / protein-group IDs are first-class keys**; never silently collapse multi-protein groups.
- **D11 is primary human evidence; D12 is validation only** (no double counting).
- CSF-specific controls: **blood contamination flags**, **secretion/localisation annotation** (including secretory autophagy cargoes), **tissue plausibility**.
- **Peptide/assay feasibility** promoted to ranking-relevant step.
- **Module (co-abundance) validation** on Astral data as independent biological check.

---

## 1. Goals and Deliverables

### Primary deliverables
1. **Ranked candidate table** (`candidates_ranked.tsv`) with per-protein columns:
   - Mouse CSF evidence (breadth, detection fraction, intensity stats, confidence tier)
   - Human CSF evidence (D11 discovery; D12 validation)
   - Autophagy/lysosome membership (R1 categories + weight)
   - EV support (annotation only; rank percentile where available)
   - Contamination flags, secretion route, tissue plausibility
   - Peptide/assay feasibility + cross-species peptide conservation
   - Astral module membership
   - Composite evidence score with all component sub-scores

2. **Core panel shortlist** (`core_panel_shortlist.tsv`) — configurable size (default top 30–80), chosen by score + feasibility tier.

3. **QC and validation package:**
   - R1 autophagy list validation against GO/Reactome (gaps and discrepancies)
   - Sensitivity analyses across parameter grid
   - Chance expectation simulation (null convergence rate)
   - GO/Reactome enrichment of convergent-but-not-autophagy proteins
   - Module validation summary

4. **Reproducibility artifacts:** config snapshot, run logs, environment spec, intermediate standardised matrices, methods text draft (for manuscript).

---

## 2. Project Layout

```
CSF_panel_project/
├── DataProc/
│   ├── Archived/                          # Previous plan versions
│   ├── DerivedData/
│   │   ├── standardised/                  # Per-dataset standardised CSVs (Step 1)
│   │   ├── orthology_cache.tsv            # Cached orthology mappings (Step 2)
│   │   ├── uniprot_cache/                 # Cached UniProt sequences (Step 5)
│   │   ├── evidence_scores/               # Scored protein lists (Step 3)
│   │   ├── autophagy_panel/               # Filtered + annotated panel (Step 4)
│   │   ├── peptide_feasibility/           # Digest and conservation results (Step 5)
│   │   └── modules/                       # WGCNA/correlation module results (Step 6)
│   ├── Log/                               # Run logs + config snapshots
│   ├── Outputs/                           # Final deliverables only
│   │   ├── candidates_ranked.tsv
│   │   ├── core_panel_shortlist.tsv
│   │   ├── overview_matrix_filled.xlsx
│   │   ├── module_validation_summary.tsv
│   │   ├── sensitivity_analysis_summary.tsv
│   │   ├── ad_model_crosscheck.tsv
│   │   ├── ev_gate_shortlist.tsv
│   │   ├── ev_gate_comparison.txt
│   │   ├── master_pipeline_list.tsv
│   │   ├── supplementary_r1b_edison_cosmos_rescored.tsv
│   │   ├── supplementary_r1b_edison_cosmos_core_qualified.tsv
│   │   ├── supplementary_r1b_edison_cosmos_not_detected.tsv
│   │   ├── supplementary_r1b_edison_cosmos_combined_ranking.tsv
│   │   ├── supplementary_r1b_edison_cosmos_report.txt
│   │   ├── figures/
│   │   └── methods_draft.md
│   ├── QC/                                # All QC artifacts
│   │   ├── contaminant_report.tsv
│   │   ├── R1_validation_vs_GO_Reactome.tsv
│   │   ├── sample_missingness_flags.tsv
│   │   ├── null_simulation_results.tsv
│   │   ├── convergent_non_autophagy_enrichment.tsv
│   │   └── ad_model_crosscheck_summary.txt
│   └── Scripts/
│       ├── config.yaml
│       ├── environment.yml
│       ├── 00_setup.py
│       ├── 01_extract_and_qc.py
│       ├── 02_orthology_mapping.py
│       ├── 03_evidence_scoring.py
│       ├── 04_autophagy_filter.py
│       ├── 05_peptide_feasibility.py
│       ├── 06_module_validation.py
│       ├── 07_validate_and_report.py
│       ├── 08_ad_model_crosscheck.py
│       ├── 09_ev_gate_analysis.py
│       ├── 10_supplementary_gene_scoring.py
│       └── utils/
│           ├── parsers.py
│           ├── orthology.py
│           ├── scoring.py
│           ├── qc.py
│           └── visualization.py
│
└── raw/                                   # Read-only source data
    ├── Full list with autophagy proteins curated.xlsx
    ├── Pooled datasets blasting overview.xlsx
    ├── Project Astral/
    │   ├── Astral protome_discovery and repliction_DACboth013_49D_1_7_3_Sc01.tsv
    │   ├── Astral protome_discovery and repliction_DACboth013_49D_1_7_3_Sc01.tsv.xlsx
    │   └── Astral proteome_repliction only_DACRP014_49D_1_7_2_Sc03.tsv
    ├── Project PXD032782/
    │   ├── APPPS1_BI_2_1Pep.xlsx
    │   └── APPPS1_BI_2_2Pep.xlsx
    ├── Project PXD048864/
    │   └── report.pg_matrix_CSF.tsv
    ├── Project PXD052590/
    │   └── report.pg_matrix.tsv
    ├── Project PXD064570/                 # Empty placeholder
    ├── Project PXD065438/
    │   ├── Output_Brain-mk_Trem2/
    │   │   └── AGCH2501_brain.pg_matrix.tsv
    │   └── Output_CSF-mk_Trem2/
    │       └── CSF-mk_Trem2.pg_matrix.tsv
    ├── The list of all identified proteins from lable-free MS of mouse CSF(Per).xlsx
    ├── WT vs KO EV protein list from highest to lowest abundance.xlsx
    ├── Table S3 Alterations in proteins identifed in all samples of AppNL-G-F NL-G-F vs Appwt wt.xlsx
    └── Misc/
        ├── CSF_Panel_Pathway_Additions.xlsx
        ├── TDP43_additional_biomarkers.xlsx
        └── novel_candidates_INCLUDE_56_extra_from_Edison_Cosmos_run.xlsx
```

---

## 3. Inputs and Dataset Registry

### 3.1 Reference inputs

| ID | File | Contents | Use |
|----|------|----------|-----|
| R1 | `raw/Full list with autophagy proteins curated.xlsx` | 604 unique human gene symbols, 6 functional categories. Use clean deduplicated sheet ("Sheet1"). | Pathway filter + scoring weight. **Must be validated against GO/Reactome before use.** |
| R2 | `raw/Pooled datasets blasting overview.xlsx` → Sheet: `Overview` | ~109 target proteins with HPA antibody IDs. | Overlap comparison only. |

### 3.2 Proteomics datasets

| ID | Species | Biofluid | Proteins | Samples | Platform | Software | Source file | Role |
|----|---------|----------|----------|---------|----------|----------|-------------|------|
| D1 | Mouse | CSF | 1,075 | 12 | Orbitrap LFQ | MaxQuant | `raw/The list of all identified proteins from lable-free MS of mouse CSF(Per).xlsx` | Mouse CSF evidence (weak: presence-only from label-free) |
| D2 | Mouse | CSF | 1,242 | 43 | Orbitrap LFQ | MaxQuant | `raw/Pooled datasets blasting overview.xlsx` → Sheet: `mouse_csf_PXD053698` | Mouse CSF evidence |
| D3 | Mouse | CSF | 1,397 | 10 | Orbitrap LFQ | MaxQuant | `raw/Pooled datasets blasting overview.xlsx` → Sheet: `mouse_csf_PXD053568_5xFAD` | Mouse CSF evidence |
| D4 | Mouse | CSF | 3,584 | 32 | timsTOF DIA | DIA-NN | `raw/Project PXD032782/APPPS1_BI_2_2Pep.xlsx` (use 2-peptide version) | Mouse CSF evidence |
| D5 | Mouse | CSF | 2,171 | 15 | timsTOF DIA | DIA-NN | `raw/Project PXD065438/Output_CSF-mk_Trem2/CSF-mk_Trem2.pg_matrix.tsv` | Mouse CSF evidence |
| D6 | Mouse | Brain | 7,645 | 15 | Orbitrap DIA | DIA-NN | `raw/Project PXD065438/Output_Brain-mk_Trem2/AGCH2501_brain.pg_matrix.tsv` | **Plausibility layer** (CNS origin confirmation) |
| D7 | Mouse | CSF | 1,555 | 20 | Orbitrap DIA | DIA-NN | `raw/Project PXD052590/report.pg_matrix.tsv` | Mouse CSF evidence |
| D8 | Mouse | CSF | 5,961 | 49 | timsTOF DIA | DIA-NN | `raw/Project PXD048864/report.pg_matrix_CSF.tsv` | Mouse CSF evidence |
| D9 | Mouse | ISF | 789 | 47 | timsTOF DIA | DIA-NN | `raw/Pooled datasets blasting overview.xlsx` → Sheet: `mouse_isf_PDX048864` | **Plausibility layer** (ISF orthogonal support) |
| D10 | Human | CSF | 990 | 8 | Orbitrap LFQ | MaxQuant | `raw/Pooled datasets blasting overview.xlsx` → Sheet: `human_csf_PXD053698` | Human CSF supporting evidence |
| D11 | Human | CSF | 3,232 | 2,720 | Astral DIA | DIA-NN | `raw/Project Astral/Astral protome_discovery and repliction_*.tsv` | **Primary human CSF evidence** |
| D12 | Human | CSF | 2,672 | 643 | Astral DIA | DIA-NN | `raw/Project Astral/Astral proteome_repliction only_*.tsv` | **Validation only** (subset of D11; never count as independent) |
| EV | Human (SH-SY5Y) | EVs | 6,097 | Pooled | — | — | `raw/WT vs KO EV protein list from highest to lowest abundance.xlsx` | **Annotation only** (not a gate) |

**D2, D3, D9, D10** are read from the pooled workbook (no standalone raw matrix on disk). All others prefer standalone raw files.

**AD model reference (Step 8 only):**

| ID | File | Contents | Use |
|----|------|----------|-----|
| Table S3 | `raw/Table S3 Alterations in proteins identifed in all samples of AppNL-G-F NL-G-F vs Appwt wt.xlsx` | CSF proteomics from APP knock-in mice (3 sheets: AppNL-G-F vs WT, AppNL-F vs WT, AppNL-F vs AppNL-G-F). ~250 proteins with log2 LFQ intensities, fold-change ratios, and p-values per comparison. | Cross-check annotation — identifies which pipeline candidates show significant changes in an AD mouse model. |

### 3.3 Key data characteristics

**MaxQuant datasets (D1–D3, D10):**
- Identifiers: `Protein IDs` (sp|ACC|NAME format), `Fasta headers` with `GN=SYMBOL`
- Quantification: raw LFQ intensities, not log-transformed
- D1 requires FASTA header parsing (`GN=(\w+)` regex); D2/D3/D10 have `gene_symbol` columns
- Sparsity: 32–56% non-null (D1); D2/D3 use zeros instead of NaN (treat zeros as missing)

**DIA-NN datasets (D4–D9, D11, D12):**
- Identifiers: `Protein.Group` + `Genes` (semicolon-delimited; may be multi-protein groups)
- Quantification: raw intensities (D4–D9); log2-transformed (D11/D12)
- Completeness: D6 ~97%; D5 79–91%; D11 65–85%; D8 18–76%; D7 21–73%
- D12 includes precomputed stats: `abundance_median_Astral`, `data_completeness_Astral`

**EV dataset:**
- Rank-ordered gene symbols only (no intensities). WT and KO columns.
- Pre-computed sheet: 329 proteins already intersected with R1.

---

## 4. Standardisation Rules

### 4.1 Primary keys and protein-group ambiguity

- Retain **protein group / accession IDs** as primary keys throughout:
  - DIA-NN: `Protein.Group` → derive `group_id`; retain `group_members` list
  - MaxQuant: `Protein IDs` and `Majority protein IDs`
- Every row carries:
  - `is_ambiguous_group` (True if multiple accessions AND multiple distinct gene symbols)
  - `group_members` (full list of accessions)
- **Never silently take the first gene from a semicolon list.** Ambiguous groups are either expanded (all members with flag) or excluded from core panel and retained in a secondary list. Configurable in `config.yaml`.

### 4.2 Contaminant and plasma-protein handling

Applied in Step 1 before any downstream analysis:
- Remove MaxQuant decoys (`REV__`) and common contaminants (`CON__`)
- Flag known high-abundance plasma proteins as `likely_plasma_derived`:
  - Start with: albumin (ALB), haemoglobins (HBA1/HBA2/HBB), immunoglobulins (IGH*/IGL*/IGK*), fibrinogen (FGA/FGB/FGG), complement (C3/C4A/C4B/C5–C9), apolipoproteins (APOA1/APOB/APOC*/APOE), transferrin (TF), alpha-2-macroglobulin (A2M), haptoglobin (HP), keratins (KRT*)
  - Optionally extend with external plasma reference list (see §10)
- `likely_plasma_derived` proteins: excluded from core panel, retained in supplementary output, applied as scoring penalty

### 4.3 Sample-level QC

Per dataset, flag and log:
- Samples with >80% missingness (potential failed runs)
- D9 (ISF): 2 completely empty sample columns — exclude these
- Output sample QC report to `QC/sample_missingness_flags.tsv`

### 4.4 Dataset-aware detectability tiers

For each protein in each dataset, compute: `n_samples`, `n_detected`, `detected_fraction`, median intensity, dataset-specific intensity floor (10th percentile of non-zero values).

**Tier definitions (configurable in `config.yaml`):**

| Tier | Small datasets (n ≤ 20) | Large datasets (n > 20) | Use |
|------|--------------------------|-------------------------|-----|
| **A (robust)** | n_detected ≥ 2 AND above intensity floor | detected_fraction ≥ 0.10 AND above intensity floor | Core panel eligibility |
| **B (moderate)** | n_detected ≥ 1 AND above intensity floor | detected_fraction ≥ 0.02 OR above intensity floor | Extended panel |
| **C (exploratory)** | Any detection | Any detection | Sensitivity analyses only |

Always report results under all three tiers.

---

## 5. Pipeline Steps

### Step 0: Configuration and setup (`00_setup.py`)

- Parse `config.yaml`; validate all input paths exist
- Create output directory structure under `DerivedData/`, `QC/`, `Outputs/`, `Log/`
- Log environment info (package versions, timestamps, config snapshot)

### Step 1: Extract, QC, and standardise (`01_extract_and_qc.py`)

**Per dataset:**
1. Parse identifiers (MaxQuant or DIA-NN format, as specified in §3.3)
2. Apply contaminant removal and plasma flagging (§4.2)
3. Apply sample-level QC (§4.3)
4. Compute detection stats and assign detectability tiers (§4.4)
5. Apply protein-group ambiguity rules (§4.1)

**Output per dataset** to `DerivedData/standardised/`:
CSV with columns: `group_id`, `uniprot_ids`, `gene_symbol`, `group_members`, `is_ambiguous_group`, `species`, `n_detected`, `pct_detected`, `median_intensity`, `above_intensity_floor`, `detectability_tier`, `likely_plasma_derived`, `dataset_id`

### Step 2: Orthology mapping (`02_orthology_mapping.py`)

1. **Pre-curate R1**: Run orthology mapping on all 604 autophagy genes first. Output the ~20–40 genes with non-trivial 1:many or many:1 mappings (ATG8/LC3/GABARAP family, ATG4 family, TRIM proteins, etc.) to `QC/` for manual review.
2. Map all mouse gene symbols → human orthologs using `gprofiler-official` (g:Orth, Ensembl-backed)
3. Cross-validate ambiguous cases with `mygene`
4. Policy: 1:1 → accept; 1:many → expand all with `orthology_ambiguous=True`; many:1 → flag
5. Cache all mappings to `DerivedData/orthology_cache.tsv` (avoid redundant API calls)
6. Convert all protein lists to human gene symbol namespace

### Step 3: Evidence scoring (`03_evidence_scoring.py`)

**Compute per-protein evidence score (continuous, not Boolean):**

```
score = w1 × mouse_CSF_evidence
      + w2 × human_CSF_evidence
      + w3 × EV_evidence
      + w4 × brain_plausibility
      + w5 × autophagy_membership
      + penalties
```

**Component definitions:**

- **mouse_CSF_evidence** (w1=0.25): `(datasets_detected / total_mouse_CSF_datasets) × mean(pct_detected across detected datasets)`. Uses D1–D5, D7, D8. D1 contributes presence-only (binary). Range 0–1.
- **human_CSF_evidence** (w2=0.30): Detection tier in D11 (Tier A=1.0, B=0.7, C=0.3, absent=0). D10 adds small bonus (+0.1 if detected). D12 used only as validation flag, not scored. Range 0–1.1.
- **EV_evidence** (w3=0.05): `ev_present × (1 - ev_rank_percentile)` where lower rank = higher abundance = higher score. **Annotation only — small weight.** Range 0–1.
- **brain_plausibility** (w4=0.10): Detected in D6 brain lysate = 1.0; absent = 0. Range 0–1.
- **autophagy_membership** (w5=0.25): In R1 = base 1.0, with category modifiers: core machinery/lysosome = 1.0; mitophagy/docking = 0.85; upstream/regulators = 0.7. Not in R1 = 0. Range 0–1.
- **penalties**: −0.1 if `is_ambiguous_group` unresolved; −0.1 if `orthology_ambiguous` unresolved; −0.2 if `likely_plasma_derived`.

Default weights configurable in `config.yaml`. Score range approximately 0–1.

**Additionally output:**
- Full provenance matrix: per-protein × per-dataset detection flag + tier
- Detection breadth: number of datasets, total samples detected
- ISF support (D9): binary annotation column
- UpSet plot: mouse CSF (pooled), human CSF, EV, brain, autophagy list

### Step 4: Autophagy filtering and biological annotation (`04_autophagy_filter.py`)

**R1 validation (run once, output to QC/):**
- Cross-reference R1 against GO:0006914 descendants and Reactome R-HSA-9612973
- Report: genes in R1 but not GO/Reactome (novel or miscurated?); genes in GO/Reactome autophagy but not R1 (gaps — especially selective autophagy receptors, ER-phagy, aggrephagy, TFEB/TFE3 regulators)
- Output: `QC/R1_validation_vs_GO_Reactome.tsv`

**Biological annotations per candidate:**
- Functional category from R1 "Full list" sheet
- Subcellular localisation: lysosome/autophagosome/endosome/cytosol (UniProt keywords or GO terms)
- Secretion route:
  - Signal peptide (classical secretion)
  - Non-classical secretion (SecretomeP prediction or known unconventional)
  - Transmembrane
  - Known EV association
  - **Known secretory autophagy cargo** (SEC22B/TRIM16-dependent pathway, IL-1β-like unconventional secretion, annexins, ferritin) — flag specifically as `secretory_autophagy_cargo`
- Brain expression plausibility: D6 detection + HPA brain tissue expression if available via API

**Enrichment of convergent-but-not-autophagy proteins:**
- Take all proteins with mouse + human CSF evidence (score components) but NOT in R1
- Run GO/Reactome enrichment analysis
- Report enriched terms — may reveal lysosome/endosome/UPS pathway components missing from R1
- Output: `QC/convergent_non_autophagy_enrichment.tsv`

**Define core panel candidates (configurable):**
- Human CSF Tier A or B in D11
- Mouse CSF Tier A or B in ≥1 mouse dataset
- In R1 autophagy/lysosome list
- Not `likely_plasma_derived`
- EV is **never** required

### Step 5: Peptide feasibility and cross-species conservation (`05_peptide_feasibility.py`)

**Run for all core panel candidates (at minimum top N by score):**

1. Fetch UniProt sequences (human + mouse ortholog). Cache in `DerivedData/uniprot_cache/`.
2. In silico tryptic digest: missed cleavages 0–2, peptide length 7–25 aa
3. Per candidate compute:
   - `n_proteotypic_peptides_human` (unique to this protein in human proteome)
   - `n_proteotypic_peptides_mouse`
   - `n_conserved_proteotypic_peptides` (identical or near-identical sequence in both species)
   - Flag: signal peptides, transmembrane domains, known glycosylation/PTM sites overlapping candidate peptides
4. Assign **feasibility tier**:
   - **Tier I**: ≥2 conserved proteotypic peptides (ideal for translational PRM/MRM)
   - **Tier II**: Human-only proteotypic peptides available (human assay feasible; separate mouse assay needed)
   - **Tier III**: No good proteotypic peptides (not assay-viable without alternative strategy)
5. Apply feasibility as scoring modifier: Tier III candidates flagged in final output; optionally penalised in score.

### Step 6: Module / co-abundance validation (`06_module_validation.py`)

**Purpose:** Verify that the panel reflects coherent CSF pathway biology, not just a curated gene list applied to detectable proteins.

On D11 (N=2,720):
1. Build log2 intensity matrix (proteins × samples). Handle missingness (filter to proteins with ≥30% completeness, or impute with care — document choice).
2. Derive co-abundance modules (WGCNA or correlation-based clustering).
3. Test R1 gene enrichment per module (Fisher's exact or hypergeometric).
4. Output:
   - Module assignment per protein
   - Enriched "autophagy/lysosome modules" with member lists
   - **Module-derived candidates**: hub proteins in enriched modules that are NOT in R1 — these are data-driven additions to review
   - Flag panel candidates that do NOT cluster with any autophagy-enriched module (they may be individually annotated autophagy genes but don't co-vary with the pathway in CSF)

### Step 7: Validation, sensitivity, and reporting (`07_validate_and_report.py`)

**Comparisons:**
- New panel vs R2 (109-target antibody panel): confirmed / new / missed
- New panel vs previous 329-protein intersection: reproducibility check
- Fill Overview template: per-target, per-dataset detection matrix → `Outputs/overview_matrix_filled.xlsx`

**Sensitivity analyses (output to `QC/` and `Outputs/sensitivity_analysis_summary.tsv`):**

Run pipeline under parameter grid:
- Detectability tier: A-only / A+B / A+B+C
- Protein-group handling: expand ambiguous / exclude ambiguous
- D4 peptide filter: 1-peptide vs 2-peptide input
- Evidence score weights: default / equal weights / CSF-heavy (w1+w2 = 0.7)
- Report: which candidates appear under all/most parameter combinations (stable core)
- Compute **`inclusion_count`** per protein: number of 9 configs (3 tier × 3 weight) where the protein ranks in the top 80. Added to `candidates_ranked.tsv` and `core_panel_shortlist.tsv` as a robustness metric (range 0–9).

**Chance expectation simulation:**
- Randomly sample 1,000 gene sets of size 604 from the human proteome (~20,000 genes)
- For each random set, compute how many genes pass the same mouse CSF + human CSF convergence criteria
- Report observed autophagy convergence count vs null distribution (p-value, fold enrichment)
- Output: `QC/null_simulation_results.tsv`

**Figures (to `Outputs/figures/`):**
- UpSet plot: mouse CSF / human CSF / EV / brain / R1 overlaps
- Detection heatmap: panel candidates × datasets
- Evidence score distribution: autophagy panel vs non-autophagy convergent
- Sensitivity stability: candidate rank consistency across parameter grid
- Module enrichment: autophagy-enriched modules visualised
- Comparison Venn: new panel vs R2 vs previous 329-protein list

**Methods text draft** → `Outputs/methods_draft.md`

**Data dictionary** → `Outputs/data_dictionary.csv`: column-level documentation (file, column, type, description) for all output files.

### Step 8: AD model cross-check (`08_ad_model_crosscheck.py`)

**Purpose:** Cross-reference pipeline candidates against CSF proteomics from APP knock-in Alzheimer's disease mouse models (Table S3) to identify disease-relevant changes among the panel proteins. This is an annotation/cross-check step — it does not modify composite scores.

**Input:** Table S3 from the AppNL-G-F study — CSF proteomics with fold changes and p-values for ~250 proteins across three genotype comparisons:
1. AppNL-G-F/NL-G-F vs Appwt/wt (full knock-in with Swedish + Iberian + Arctic mutations vs wild-type)
2. AppNL-F/NL-F vs Appwt/wt (Swedish + Iberian mutations vs WT)
3. AppNL-F/NL-F vs AppNL-G-F/NL-G-F (partial vs full knock-in)

**Processing:**
1. Parse mouse gene symbols from UniProt FASTA-header entries (`GN=` tags via regex)
2. Map to human orthologs via the Step 2 orthology cache (fallback: uppercase match)
3. Cross-reference against `candidates_ranked.tsv` (all 9,561 scored proteins)
4. Annotate each overlapping protein with ratio, p-value, and direction per comparison
5. Classify direction: `up` (ratio>1, p<0.05), `down` (ratio<1, p<0.05), `ns` (p≥0.05)
6. Compute summary columns: `n_comparisons_significant` (0–3), `consistent_direction` (bool)

**Output:**
- `Outputs/ad_model_crosscheck.tsv` — annotated overlap table with per-comparison fold changes, p-values, directions, plus summary columns
- `QC/ad_model_crosscheck_summary.txt` — human-readable report of overlaps, significant hits in the top 80 shortlist, and all overlapping proteins

### Step 9: EV hard-gate sensitivity analysis (`09_ev_gate_analysis.py`)

**Purpose:** Test the effect of making EV support (`ev_present`) a hard pass/fail requirement for core panel membership, instead of just a scoring weight. This is a read-only analysis — it loads existing pipeline outputs and re-applies gating logic in-memory. No existing output files are modified.

**Processing:**
1. Load `candidates_ranked.tsv` and `core_panel_shortlist.tsv` (existing outputs)
2. Apply all 4 current gates plus a 5th EV gate (`ev_present=True`)
3. Build a new top-80 shortlist from the surviving candidates
4. Compare against the current shortlist: proteins dropped, gained, and rank shifts
5. Cross-reference dropped proteins against AD model data (if available from Step 8)

**Output:**
- `Outputs/ev_gate_shortlist.tsv` — top 80 under EV hard-gate (same column structure as `core_panel_shortlist.tsv`)
- `Outputs/ev_gate_comparison.txt` — human-readable diff report: gate cascade, dropped/gained proteins, rank shifts, AD model overlap for dropped proteins

### Step 10: Supplementary gene scoring — R1b (`10_supplementary_gene_scoring.py`)

**Purpose:** Score supplementary curated gene lists (designated "R1b") against the existing pipeline evidence. As the curated autophagy/lysosome gene list grows over time, new gene sets can be plugged in via `config.yaml` without re-running the core pipeline. The module recalculates autophagy membership scores for genes that were already detected in the pipeline's datasets but were not in the original R1 reference list.

**Design:** Reusable — parameterised via `config.yaml` `supplementary_gene_lists` section. Each entry specifies an Excel file with gene symbols, functional categories, and database sources. Multiple supplementary lists can be processed in a single run.

**First list:** `R1b_edison_cosmos` — 56 novel autophagy genes + 9 HADb core proteins from a systematic Edison/Cosmos database search (`raw/Misc/novel_candidates_INCLUDE_56_extra_from_Edison_Cosmos_run.xlsx`). None overlap with the original R1 list (599 genes).

**Processing:**
1. Load supplementary gene list from Excel; expand slash-separated gene symbols
2. Build category map: HADb Core sheet provides functional categories for 9 genes (Mitophagy, CMA, fusion, etc.); remaining genes inferred from Database Sources keywords (GO/KEGG/Reactome terms); fallback to "default" tier
3. Skip genes already in R1 (TECPR1, EPG5 — already scored in original run)
4. Load `evidence_scores.tsv` (all 9,561 proteins with pre-computed score components)
5. For found genes: recalculate `score_autophagy` = 1.0 × category_modifier (was 0), recompute `composite_score` with original weights (sum=1.0), clamp [0,1]
6. Apply all 4 hard gates (D11 tier A/B, mouse tier A/B, in_r1b=True, not plasma)
7. Build combined ranking: merge R1 core (122) + R1b qualified, sort by composite_score

**Output:**
- `Outputs/supplementary_r1b_edison_cosmos_rescored.tsv` — all found genes with old vs new scores
- `Outputs/supplementary_r1b_edison_cosmos_core_qualified.tsv` — genes passing all 4 gates
- `Outputs/supplementary_r1b_edison_cosmos_not_detected.tsv` — genes absent from all datasets
- `Outputs/supplementary_r1b_edison_cosmos_combined_ranking.tsv` — R1 + R1b merged ranking
- `Outputs/supplementary_r1b_edison_cosmos_report.txt` — human-readable summary

---

## 6. Output Contract

**Every row in every intermediate and final output must retain:**
- `dataset_id` (where applicable)
- `species`
- Original identifiers (`uniprot_ids`, `group_id`, `group_members`)
- Mapped `human_gene_symbol`
- `is_ambiguous_group` flag
- `orthology_ambiguous` flag
- `detectability_tier` (per dataset)
- All evidence score components (not just the composite)

No information is silently discarded. Filtering decisions are recorded as flag columns, not row deletions, until the final shortlist step.

---

## 7. External Data Access

Pipeline runs offline with `raw/`. Optional additions to strengthen specific evidence layers:

### 7.1 CSF-derived EV proteomics (recommended)
If EV evidence needs upgrading beyond the SH-SY5Y cell-line data:
- Create `raw/External_CSF_EV/` and add published CSF-EV proteomics tables (e.g., Chiasserini et al., Muraoka et al.)
- Record provenance in a `README.md` in that directory
- Use as improved EV annotation layer (still not a gate)

### 7.2 Plasma reference lists
- Create `raw/External_PlasmaRefs/` with a TSV of plasma-enriched proteins (gene symbols + UniProt accessions)
- Improves contamination flagging beyond the built-in curated list
- Record provenance

### 7.3 Missing raw matrices
- `raw/Project PXD064570/` is currently empty
- If needed: download processed quantitative matrices (not vendor RAW files) from PRIDE/ProteomeXchange
- Save to corresponding `raw/Project <ACCESSION>/` directory

---

## 8. Python Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| `pandas` | ≥2.0 | Data manipulation |
| `openpyxl` | ≥3.1 | Excel I/O |
| `gprofiler-official` | ≥1.0 | Orthology mapping (g:Orth) |
| `mygene` | ≥3.2 | Gene annotation, cross-validation |
| `pyteomics` | ≥4.6 | In silico tryptic digest |
| `biopython` | ≥1.83 | Sequence retrieval, alignment |
| `gseapy` | ≥1.0 | GO/Reactome enrichment analysis |
| `matplotlib` | ≥3.7 | Plotting |
| `matplotlib-venn` | ≥0.11 | Venn diagrams |
| `upsetplot` | ≥0.9 | Multi-set UpSet plots |
| `seaborn` | ≥0.13 | Heatmaps |
| `pyyaml` | ≥6.0 | Config parsing |
| `requests` | — | UniProt/HPA REST API |
| `scipy` | ≥1.11 | Statistical tests, clustering |

Optional for module analysis: `rpy2` + R `WGCNA`, or pure-Python correlation clustering.

---

## 9. Implementation Notes (for Claude Code)

- **Single config file**: `DataProc/Scripts/config.yaml` listing all dataset IDs, file paths, format specs (MaxQuant vs DIA-NN), key column mappings, all tuneable parameters (weights, thresholds, tier definitions).
- **Deterministic transforms**: no randomness without seed; always write intermediates; no hidden state.
- **Stepwise execution**: each script reads prior step's output from `DerivedData/` and writes to the next. Scripts can be run independently after Step 0.
- **Logging**: every script logs start/end, input/output paths, row counts, and parameter values to `Log/`.
- **Fail loudly**: if an input file is missing or a parsing step produces zero rows, raise an error — do not silently continue with empty data.
- **Output contract** (§6): enforced in every write operation.
