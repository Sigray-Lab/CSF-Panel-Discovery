# Methods Draft: CSF Autophagy/Lysosome Panel Discovery

## Proteomics Data Processing

We analyzed 12 proteomics datasets spanning mouse CSF (7 datasets),
mouse brain tissue (1 dataset), mouse ISF (1 dataset), and human CSF (3 datasets),
comprising a total of 9561 unique proteins after orthology mapping to human
gene symbols.

Raw proteomics data from MaxQuant (label-free quantification) and DIA-NN
(data-independent acquisition) pipelines were standardised using a uniform
processing pipeline. Contaminants (decoy sequences, common laboratory
contaminants) and high-abundance plasma proteins were flagged and excluded
from the core panel. Sample-level quality control flagged samples with
>80% missing values.

## Detectability Tiers

Proteins were assigned dataset-aware detectability tiers (A/B/C) based on
detection frequency and intensity relative to a dataset-specific intensity
floor (10th percentile of non-zero values). For small datasets (n <= 20
samples), Tier A required detection in >= 2 samples above the intensity
floor. For large datasets (n > 20), Tier A required >= 10% detection
fraction above floor.

## Orthology Mapping

Mouse gene symbols were mapped to human orthologs using g:Profiler g:Orth
(Ensembl-backed) with cross-validation via MyGene.info for ambiguous cases.
One-to-one orthologs were accepted directly; one-to-many mappings were
expanded with an ambiguity flag.

## Evidence Scoring

A continuous composite evidence score was computed for each protein:

Score = 0.25 x mouse_CSF + 0.3 x human_CSF
      + 0.1 x EV + 0.1 x brain
      + 0.25 x autophagy + penalties

Mouse CSF evidence was computed as (proportion of datasets with detection)
x (mean detection fraction across detected datasets). Human CSF evidence
was based on the Astral discovery cohort (D11, N=2,720) detectability tier.
EV evidence served as supportive annotation only (weight = 0.1).
Brain tissue plausibility was assessed from mouse brain lysate proteomics.
Autophagy/lysosome membership was scored based on a curated reference list
of 434 genes across 6 functional categories.

## Core Panel Selection

Core panel candidates were required to meet: (1) Human CSF detectability
tier A or B in the Astral cohort, (2) Mouse CSF tier A or B in at least
one dataset, (3) Membership in the curated autophagy/lysosome reference
list, and (4) Not flagged as likely plasma-derived.

## Sensitivity Analysis

Pipeline robustness was assessed across a parameter grid of 18 combinations
varying detectability tier stringency (A-only, A+B, A+B+C), protein group
handling (expand vs exclude ambiguous), and evidence score weight sets
(default, equal, CSF-heavy).

## Chance Expectation

To assess the significance of autophagy gene convergence across species,
we performed 1,000 random simulations sampling 604 genes from the human
proteome and computing how many passed the same mouse + human CSF
convergence criteria.

## Software

The pipeline was implemented in Python (>= 3.10) using pandas, scipy,
gprofiler-official, mygene, pyteomics, biopython, gseapy, matplotlib,
seaborn, and upsetplot. All analyses were deterministic (seed = 42) with
intermediate results saved at each step.
