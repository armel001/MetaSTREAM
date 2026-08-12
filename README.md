<div align="center">

# MetagenAMR

**Shotgun Nanopore metagenomics pipeline for antimicrobial resistance surveillance**

[![Snakemake](https://img.shields.io/badge/Snakemake-≥7.0-brightgreen?style=flat-square&logo=python)](https://snakemake.readthedocs.io)
[![Conda](https://img.shields.io/badge/Conda-supported-44A833?style=flat-square&logo=anaconda)](https://docs.conda.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue?style=flat-square)](LICENSE)
[![Nanopore](https://img.shields.io/badge/Oxford_Nanopore-R10.4-0084C8?style=flat-square)](https://nanoporetech.com)
[![Status](https://img.shields.io/badge/Status-Active-success?style=flat-square)]()

*Characterization of resistomes, mobilomes, and microbial communities from Oxford Nanopore long-read sequencing data*

</div>

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Workflow](#pipeline-workflow)
- [Requirements](#requirements)
- [Installation](#installation)
- [Configuration](#configuration)
- [Quick Test](#quick-test)
- [Usage](#usage)
- [Outputs](#outputs)
- [Related Projects](#related-projects)
- [Citation](#citation)
- [License](#license)

---

## Overview

MetagenAMR is a reproducible, end-to-end Snakemake pipeline integrating quality control, taxonomic profiling, resistome and mobilome analysis from Oxford Nanopore shotgun metagenomics data. It produces publication-ready figures and structured output matrices suitable for downstream statistical analysis.

**Key features:**
- Long-read aware normalization (copies/Gb using per-sample mean read length)
- Dual ARG detection strategy — assembly-based (RGI main) and read-based (RGI BWT)
- Contig-level ARG × MGE co-localization analysis
- Automated HTML reports for all analytical modules
- Fully reproducible via Conda environments per rule

---

## Pipeline Workflow

```
Raw reads ─▶ QC & Assembly ─▶ Taxonomy ─▶ Resistome (RGI) ─┐
                                                             ├─▶ Reports & Matrices
                                          Mobilome (MGE) ────┘
```

---

## Requirements

| Requirement | Version | Notes |
|---|---|---|
| [Snakemake](https://snakemake.readthedocs.io) | ≥ 7.0 | Workflow manager |
| [Conda](https://docs.conda.io) / Mamba | any | Environment management |
| Python | ≥ 3.10 | Pipeline scripts |
| R | ≥ 4.3 | Diversity & normalization |
| Disk space | ~50 GB/sample | Assembly + intermediate files |

```bash
conda install -n base -c conda-forge snakemake mamba
```

---

## Installation

```bash
git clone https://github.com/armel001/MetagenAMR.git
cd MetagenAMR
```

All tool environments are automatically built by Snakemake via `--use-conda` — no manual installation required.

**Databases.** MetagenAMR relies on external reference databases (Kraken2, CARD/RGI, mobileOG-db, PlasmidFinder, a host reference genome) that are not distributed with the repository — they must be downloaded and placed under `resources/` following the layout described in the project documentation. Large databases can be symlinked rather than copied.

---

## Configuration

Pipeline behaviour is controlled by `config/config.yaml` (samples, scientific parameters, database paths). Two optional overlays:
- `config/config.local.yaml` — machine-specific settings (never committed).
- `config/config.test.yaml` — runs the pipeline on the small bundled test dataset (see below).

```bash
snakemake --configfile config/config.yaml config/config.local.yaml -n
```

---

## Quick Test

A small bundled dataset lets you validate the installation end-to-end without real samples:

```bash
cp resources/data_test/*.fastq.gz resources/reads/
snakemake --configfile config/config.yaml config/config.test.yaml --use-conda -n
```

This is a technical smoke test (confirms the pipeline runs), not a biological validation.

---

## Usage

```bash
bash run_AMR_pipeline.sh          # complete pipeline
bash run_AMR_pipeline.sh -t <target>
bash run_AMR_pipeline.sh -n       # dry-run
```

Main targets: `taxonomy_all`, `rgi_main_all`, `rgi_bwt_all`, `mge_all`, `cooccurrence_all`, and `all` (complete pipeline). Run `snakemake --list` for the exhaustive list of rules.

**Note on Medaka polishing:** requires a GPU. If unavailable locally, polishing can be run on a remote GPU server and results copied back into `results/{sample}/medaka/` — see project documentation.

---

## Outputs

Per-sample results under `results/{sample}/` (QC, assembly, RGI, MGE hits); pooled matrices under `results/r_analysis/`; HTML reports under `results/figures/`. Full output map: see project documentation.

---

## Related Projects

- **MAG_pipeline** — companion pipeline for metagenome-assembled genome (MAG) reconstruction and quality assessment, developed for the same study. *Repository coming soon.*

---

## Citation

If you use MetagenAMR in your research, please cite:

> Gnimadi TAC et al. (2026) Manuscript in preparation.

Built on Snakemake, NanoPlot, Minimap2, fastplong, Flye, Kraken2/Bracken, RGI/CARD, DIAMOND + mobileOG-db, and PlasmidFinder — please also cite these dependencies where relevant. Full versions and citations: see project documentation.

---

## License

MIT License — Copyright (c) 2026 Thibaut Armel Chérif GNIMADI et al. See [LICENSE](LICENSE) for full terms.

---

<div align="center">
<sub>MetagenAMR · Snakemake · Conda · Nanopore R10.4</sub>
</div>
