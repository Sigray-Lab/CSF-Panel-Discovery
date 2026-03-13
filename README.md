# CSF Autophagy & Lysosome Biomarker Panel Discovery

A reproducible, config-driven proteomics pipeline that identifies and ranks **cerebrospinal fluid (CSF)-detectable markers of secretory autophagy and lysosomal function** with cross-species translation, integrating **12 mass-spectrometry datasets across ~3,600 samples**.

> **Which autophagy and lysosomal proteins can we reliably measure in CSF across both mouse and human?**
>
> We need candidates that are (1) robustly detected in human and mouse CSF, (2) biologically relevant to secretory autophagy, macro-autophagy and lysosomal pathways, and (3) directly amenable to targeted assay development with shared peptide sequences across species.

---

## Key Results

| | |
|---|---|
| **9,561** unique proteins scored | across 5 weighted evidence axes |
| **122** core panel candidates | pass all hard gates (human CSF + mouse CSF + autophagy membership + plasma exclusion) |
| **80** shortlisted proteins | ranked by composite evidence score |
| **99/122** Tier I assay-ready | shared proteotypic peptides for direct mouse-human PRM/MRM bridging |
| **21/80** module-validated | co-vary with known autophagy markers in independent co-expression analysis |

---

## Pipeline Overview

```
 12 datasets (~3,600 samples)
 ───────────────────────────────────────────────────────────────────
  Mouse CSF (7)  |  Mouse brain (1)  |  Human CSF (3)  |  EV ref (1)
 ───────────────────────────────────────────────────────────────────
         │                                     │
    01  Extract & QC ──────────────────── Extract & QC
         │                                     │
    02  Orthology mapping (g:Profiler)         │
         │       89% clean 1:1 mappings        │
         └──────────────┬──────────────────────┘
                        │
    03  Evidence scoring (5-axis weighted composite)
                        │
    04  Autophagy hard-gate filter ──► 122 core candidates
                        │
    05  Peptide feasibility (cross-species conservation)
                        │
    06  Co-abundance module validation (Astral, N=2,720)
                        │
    07  Sensitivity analysis & QC ──► 80-protein shortlist
                        │
    08  AD model cross-check (App knock-in mice)
    09  EV gate sensitivity analysis
    10  Supplementary gene list scoring
```

---

## Input Data

### Mouse datasets (9)

| ID | Biofluid | Proteins | Samples | Platform | Role |
|----|----------|----------|---------|----------|------|
| D1 | CSF | 1,075 | 12 | Orbitrap LFQ | Mouse CSF evidence |
| D2 | CSF | 1,242 | 43 | Orbitrap LFQ | Mouse CSF evidence |
| D3 | CSF | 1,397 | 10 | Orbitrap LFQ | Mouse CSF evidence |
| D4 | CSF | 3,584 | 32 | timsTOF DIA | Mouse CSF evidence |
| D5 | CSF | 2,171 | 15 | timsTOF DIA | Mouse CSF evidence |
| D6 | Brain | 7,645 | 24 | timsTOF DIA | Brain plausibility |
| D7 | CSF | 5,961 | 12 | timsTOF DIA | Mouse CSF evidence |
| D8 | CSF | 3,399 | 40 | timsTOF DIA | Mouse CSF evidence |
| D9 | ISF | 789 | 15 | Orbitrap LFQ | Orthogonal support |

### Human datasets (3 + references)

| ID | Biofluid | Proteins | Samples | Platform | Role |
|----|----------|----------|---------|----------|------|
| D10 | CSF | 990 | 49 | DIA | Independent validation |
| D11 | CSF | 3,232 | 2,720 | Astral DIA | **Primary human evidence** |
| D12 | CSF | (replication subset of D11) | | Astral DIA | Replication validation |
| R1 | -- | 604 genes | -- | -- | Curated autophagy/lysosome list |
| EV | -- | SH-SY5Y secretome | -- | -- | EV annotation (5% weight) |

---

## Methods

### Quality Control

- **Sample filtering**: Samples with >80% missing data excluded; MaxQuant decoys and lab contaminants removed.
- **Blood contamination**: 235 proteins flagged as likely plasma-derived (albumin, immunoglobulins, complement, keratins) and excluded from the core panel.
- **Detectability tiers**: Dataset-aware confidence grading scales with sample size. Tier A = robust detection; Tier B = moderate. Thresholds adapt from n=8 to n=2,720.

### Orthology Mapping

All mouse gene symbols mapped to human equivalents via g:Profiler (Ensembl gene trees):

- **89%** clean 1:1 mappings (8,589 genes)
- **5%** ambiguous 1:many (477 genes) -- flagged, kept, penalised in scoring
- **6%** no human ortholog (~604 genes) -- dropped

Paralog families (~30 core autophagy genes with non-trivial mappings, e.g. LC3/GABARAP, ATG4 variants) were identified upfront and handled explicitly. Ambiguous cases use an expand policy: all possible human matches kept but flagged.

### Evidence Scoring

Each protein receives a continuous composite score (0--1) from five weighted components:

| Component | Weight | Source |
|-----------|--------|--------|
| Human CSF evidence | 0.30 | Astral discovery cohort (D11) + validation (D10, D12) |
| Mouse CSF evidence | 0.25 | 7 mouse CSF datasets (D1--D5, D7--D8) |
| Autophagy membership | 0.25 | Curated autophagy/lysosome gene list (R1, 604 genes) |
| Brain plausibility | 0.10 | Mouse brain lysate (D6) |
| EV support | 0.05 | Human cell-line EV secretome |

No single data source determines inclusion. EV evidence contributes a small scoring bonus but is **never required**.

### Core Panel Selection (Hard Gates)

A protein must pass **all four** criteria to enter the core panel:

1. Tier A or B in human CSF (Astral, D11)
2. Tier A or B in at least one mouse CSF dataset
3. Member of the curated autophagy/lysosome reference list (R1)
4. Not flagged as likely plasma-derived

This filtering reduces 9,561 scored proteins to **122 core candidates**, ranked by composite score to yield the **top-80 shortlist**.

```
9,561 proteins scored
  └─► 6,684  detected in >=1 mouse CSF dataset
       └─► 3,026  Tier A/B in human Astral
            └─► 434  on R1 autophagy list
                 └─► 122  pass all 4 hard gates
                      └─► 80  final shortlist (by score)
```

---

## Top 20 Candidates

| Rank | Protein | Score | Feasibility | Mouse CSF |
|------|---------|-------|-------------|-----------|
| 1 | CTSD | 0.975 | Tier I | 7/7 |
| 2 | CTSB | 0.927 | Tier II | 7/7 |
| 3 | CST3 | 0.927 | Tier II | 7/7 |
| 4 | CPQ | 0.926 | Tier I | 7/7 |
| 5 | PRDX6 | 0.926 | Tier II | 7/7 |
| 6 | PARK7 | 0.925 | Tier I | 7/7 |
| 7 | HEXB | 0.922 | Tier I | 7/7 |
| 8 | GGH | 0.909 | Tier II | 7/7 |
| 9 | GAA | 0.902 | Tier I | 7/7 |
| 10 | LAMP1 | 0.896 | Tier I | 7/7 |
| 11 | ATP6AP1 | 0.887 | Tier I | 7/7 |
| 12 | NPC2 | 0.884 | Tier I | 7/7 |
| 13 | GM2A | 0.884 | Tier I | 7/7 |
| 14 | EPDR1 | 0.874 | Tier I | 7/7 |
| 15 | IDS | 0.873 | Tier I | 7/7 |
| 16 | CTSS | 0.872 | Tier I | 7/7 |
| 17 | CTSZ | 0.870 | Tier I | 7/7 |
| 18 | PPT1 | 0.869 | Tier I | 6/7 |
| 19 | MAN2B1 | 0.869 | Tier I | 7/7 |
| 20 | HEXA | 0.865 | Tier I | 7/7 |

**Feasibility**: Tier I = 2+ conserved proteotypic peptides (mouse & human); Tier II = human peptides available, separate mouse assay needed. **Mouse CSF** = datasets detected in (of 7 CSF datasets). The top 20 is dominated by lysosomal enzymes and autophagy effectors, with 19/20 detected in all 7 mouse CSF datasets.

---

## Validation

### Peptide Feasibility

Of the 122 core panel proteins:

- **99 (81%)** -- Tier I: 2+ conserved proteotypic peptides shared between mouse and human. These can be quantified with a single PRM/MRM assay using identical peptides in both species, directly bridging preclinical and clinical studies.
- **23 (19%)** -- Tier II: Human-specific peptides available; separate mouse assay design needed.
- **0** -- Tier III (unmeasurable).

### Co-abundance Module Validation

All 2,323 Astral-detected proteins were clustered by CSF co-variation across 2,720 individuals -- entirely independent of any curated gene list.

- **2 modules** significantly enriched for autophagy genes (FDR < 0.05)
- **21 of 80** shortlisted proteins fall within enriched modules
- **9 hub proteins** in enriched modules were *not* on R1 -- data-driven candidates that emerged from co-expression, not curation

This provides independent biological evidence: roughly a quarter of the panel doesn't just carry an autophagy label -- these proteins co-vary together in CSF, suggesting coherent pathway-level biology.

### AD Mouse Model Cross-check

The panel was cross-referenced against CSF proteomics from App^NL-G-F and App^NL-F knock-in mice vs wild-type (3 pairwise comparisons):

- **17/80** shortlisted proteins present in the AD model data
- **9** significantly altered in at least one comparison (p < 0.05)
- **Dominant pattern**: lysosomal enzymes decreased in CSF of amyloid models, consistent with impaired lysosomal secretion in AD

| Gene | Panel rank | Direction | p-value | Comparison |
|------|-----------|-----------|---------|------------|
| CTSD | 1 | decreased | 0.017 | NL-F vs NL-GF |
| CTSB | 2 | decreased | 0.017, 0.002 | NL-F, NL-F vs NL-GF |
| CPQ | 4 | decreased | 0.003 | NL-F vs WT |
| HEXB | 7 | decreased | 0.028 | NL-F vs WT |
| MAN2B1 | 19 | decreased | 0.046 | NL-F vs WT |

### Sensitivity Analysis

The pipeline was tested under 9 configurations varying tier stringency and scoring weights:

- **Mean Spearman rho = 0.89** across configurations (range 0.79--0.98)
- **65/80** candidates retained on average across all configurations
- **Stable core**: CTSD, CTSB, CST3, CPQ, PARK7, HEXB, GAA, and LAMP1 survived every configuration

### EV Gate Analysis

- **Up-weighting EV to 25%**: 72/80 proteins unchanged (rho = 0.88). 8 lysosomal enzymes swap out for vesicle trafficking proteins -- biologically expected.
- **Hard EV gate**: 122 drops to 102 candidates. Three of the lost proteins (CPQ, CTSS, CTSZ) are significantly altered in the AD model. An EV gate trades disease-validated hits for marginal EV-positive candidates.
- **Conclusion**: Lysosomal enzymes reach CSF via direct secretion, not EV packaging. Low EV weight (5%) is appropriate.

---

## Supplementary Gene List Scoring

Step 10 extends the pipeline to additional curated gene lists (R1b) scored against the same evidence framework. This adds **14 qualified R1b candidates** for a **136-protein combined ranking** (122 R1 + 14 R1b), allowing the panel to incorporate newly curated targets without re-running the full pipeline.

---

## Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 00 | `00_setup.py` | Validate inputs, create directory structure |
| 01 | `01_extract_and_qc.py` | Parse 12 datasets into standardised format; flag contaminants |
| 02 | `02_orthology_mapping.py` | Mouse-to-human orthology via g:Profiler (Ensembl gene trees) |
| 03 | `03_evidence_scoring.py` | Compute 5-axis weighted evidence scores for all proteins |
| 04 | `04_autophagy_filter.py` | Apply hard gates, rank candidates, generate core panel |
| 05 | `05_peptide_feasibility.py` | Cross-species peptide conservation & targeted assay feasibility |
| 06 | `06_module_validation.py` | Co-abundance module analysis on Astral data (N=2,720) |
| 07 | `07_validate_and_report.py` | Sensitivity analyses, null simulations, figures, methods draft |
| 08 | `08_ad_model_crosscheck.py` | Cross-check against App^NL-G-F AD mouse model CSF proteomics |
| 09 | `09_ev_gate_analysis.py` | EV hard-gate and reweighting sensitivity analysis |
| 10 | `10_supplementary_gene_scoring.py` | Score supplementary curated gene lists (R1b) |

All parameters are centralised in [`config.yaml`](DataProc/Scripts/config.yaml). The pipeline is fully config-driven with no hardcoded paths.

---

## Outputs

| File | Description |
|------|-------------|
| `candidates_ranked.tsv` | All 9,561 scored proteins with full evidence breakdown |
| `core_panel_shortlist.tsv` | Top 80 shortlisted proteins (pass all hard gates) |
| `master_pipeline_list.tsv` | 139 master list (122 core + 17 human-only additions) |
| `supplementary_r1b_*_combined_ranking.tsv` | 136 combined ranking (122 R1 + 14 R1b) |
| `supplementary_r1b_*_report.txt` | Supplementary scoring summary |
| `figures/` | Publication-ready figures (score distributions, sensitivity, modules) |
| `methods_draft.md` | Methods section draft |

---

## Getting Started

### Prerequisites

- [Conda](https://docs.conda.io/) or [Mamba](https://mamba.readthedocs.io/)
- Python >= 3.10
- Raw data files (see [`raw/README.md`](raw/README.md) for sourcing instructions)

### Setup & Run

```bash
# Clone the repository
git clone https://github.com/Sigray-Lab/CSF-Panel-Discovery.git
cd CSF-Panel-Discovery

# Create conda environment
conda env create -f DataProc/Scripts/environment.yml
conda activate csf_panel

# Place raw data files in raw/ (see raw/README.md)

# Run the full pipeline
cd DataProc/Scripts
for step in 00 01 02 03 04 05 06 07 08 09 10; do
  python ${step}_*.py
done
```

### Directory Structure

```
CSF_panel_project/
├── DataProc/
│   ├── Scripts/              # Pipeline code (11 steps + utilities)
│   │   ├── config.yaml       # Central configuration
│   │   ├── environment.yml   # Conda environment
│   │   ├── 00_setup.py ... 10_supplementary_gene_scoring.py
│   │   └── utils/            # Shared modules (parsers, scoring, orthology, QC, viz)
│   ├── Outputs/              # Final deliverables
│   ├── DerivedData/          # Intermediate files (regenerated by pipeline)
│   ├── QC/                   # Quality control artifacts (regenerated)
│   └── Log/                  # Timestamped execution logs (regenerated)
└── raw/                      # Source data (not included; see raw/README.md)
```

### Software Stack

pandas, scipy, gprofiler-official, mygene, pyteomics, biopython, gseapy, matplotlib, seaborn, upsetplot. All analyses are deterministic (seed = 42) with intermediate results saved at each step.

---

## License

TBD

## Citation

TBD
