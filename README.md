# tpp1-cln2-pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/release/python-3110/)
[![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A57.32-brightgreen.svg)](https://snakemake.readthedocs.io)

**Automated ACMG/AMP Variant Pathogenicity Classification Pipeline for TPP1 (CLN2 Disease)**

**Authors:** Lena Traczuk, Dawid Fleischer

---

## Overview

`tpp1-cln2-pipeline` is a fully reproducible, containerized bioinformatics pipeline for the automated classification of variants in the *TPP1* gene according to ACMG/AMP guidelines, implementing the Tavtigian Bayesian point framework. The pipeline integrates VEP v115, SpliceAI v1.3.1, gnomAD v4.1.0, and ClinVar for comprehensive variant annotation.

**Disease context:** Late infantile neuronal ceroid lipofuscinosis (CLN2 disease, Batten disease) — autosomal recessive lysosomal storage disorder.

---

## Results Summary

| Metric | Value |
|---|---|
| Total ClinVar Variants Analyzed | **607** |
| Pathogenic | 0 |
| Likely Pathogenic | **89** |
| VUS | 258 |
| Likely Benign | 15 |
| Benign | 245 |
| Clinically Actionable (P/LP) | **89 (14.7%)** |
| Concordance with Franklin | **75.8%** |
| Data Source | ClinVar GRCh38 (2026-03-21) |

---

## Pipeline Architecture

```mermaid
graph TD
    A[Raw VCF] --> B[BCFtools Normalize]
    B --> C[VEP v115 + REVEL + BayesDel]
    C --> D[SpliceAI v1.3.1]
    D --> E[gnomAD v4.1.0 AF]
    E --> F[ClinVar Annotation]
    F --> G[ACMG/AMP Classifier]
    G --> H[Results TSV]
    H --> I[Figures & Report]
    H --> J[Platform Comparison]
```

---

## Quick Start

```bash
# Clone repository
git clone https://github.com/lenax04/tpp1-cln2-pipeline.git
cd tpp1-cln2-pipeline

# Run full pipeline with conda environments
snakemake --use-conda --cores 8

# Or with Docker
docker build -t tpp1-cln2-pipeline .
docker run -v $(pwd):/pipeline tpp1-cln2-pipeline snakemake --cores 4
```

---

## Repository Structure

```
tpp1-cln2-pipeline/
├── Snakefile                          # Main workflow
├── config/config.yaml                 # Gene-specific parameters
├── rules/                             # Snakemake rule modules
│   ├── normalize.smk
│   ├── annotate_vep.smk
│   ├── annotate_spliceai.smk
│   ├── annotate_gnomad.smk
│   ├── annotate_clinvar.smk
│   ├── classify_acmg.smk
│   ├── compare_platforms.smk
│   └── report.smk
├── scripts/                           # Python analysis scripts
│   ├── classify_acmg.py               # Tavtigian Bayesian classifier
│   ├── annotate_gnomad.py
│   ├── annotate_clinvar_and_parse.py
│   ├── compare_platforms.py
│   ├── generate_figures.py
│   └── generate_report.py
├── envs/                              # Conda environments
├── tests/test_pipeline.py             # Unit tests
├── data/TPP1_raw.tsv                  # Real ClinVar variants (607)
├── results/                           # Analysis outputs
│   ├── TPP1/classification/
│   ├── TPP1/figures/
│   └── comparison/
├── paper/paper.md                     # Manuscript                 
├── Dockerfile
├── CITATION.cff
└── LICENSE
```

---

## ACMG/AMP Criteria Implemented

The classifier implements the Tavtigian (2020) Bayesian point system:

| Criterion | Evidence Level | Points | Description |
|---|---|---|---|
| PVS1 | Very Strong Pathogenic | +8 | Null variant (LoF) |
| PS3 | Strong Pathogenic | +4 | Functional studies |
| PM1 | Moderate Pathogenic | +2 | Critical domain |
| PM2 | Moderate Pathogenic | +1 | Absent from gnomAD |
| PP3 | Supporting Pathogenic | +1 to +2 | In silico (REVEL/SpliceAI) |
| BS1 | Strong Benign | -4 | AF > 1% |
| BA1 | Stand-alone Benign | -8 | AF > 5% |

---

## Citation

If you use this pipeline, please cite:

```bibtex
  author = {Traczuk, Lena and Fleischer, Dawid},
  title = {tpp1-cln2-pipeline: Automated ACMG/AMP Variant Classification for TPP1},
  year = {2026},
  url = {https://github.com/lenax04/tpp1-cln2-pipeline}
}
```
