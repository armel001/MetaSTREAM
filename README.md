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

## 📋 Table of Contents

- [Overview](#-overview)
- [Pipeline Workflow](#-pipeline-workflow)
- [Requirements](#-requirements)
- [Installation](#-installation)
- [Configuration](#-configuration)
- [Usage](#-usage)
- [Outputs](#-outputs)
- [Tools & Dependencies](#-tools--dependencies)
- [Citation](#-citation)
- [License](#-license)

---

## 🔭 Overview

MetagenAMR is a reproducible, end-to-end Snakemake pipeline integrating quality control, taxonomic profiling, resistome and mobilome analysis from Oxford Nanopore shotgun metagenomics data. It produces publication-ready figures and structured output matrices suitable for downstream statistical analysis.

**Key features:**
- Long-read aware normalization (copies/Gb using per-sample mean read length)
- Dual ARG detection strategy — assembly-based (RGI main) and read-based (RGI BWT)
- Contig-level ARG × MGE co-localization analysis
- Automated HTML reports for all analytical modules
- Fully reproducible via Conda environments per rule

---

## 🔄 Pipeline Workflow

```
 ┌─────────────────────────────────────────────────────────────────┐
 │                                                                 │
 │   Raw reads (.fastq.gz)                                         │
 │        │                                                        │
 │        ▼                                                        │
 │   ┌─────────────┐   ┌──────────────┐   ┌──────────────────┐   │
 │   │  01 · QC &  │──▶│ 02 · Taxo-  │──▶│  03 · Resistome  │   │
 │   │  Assembly   │   │    nomy      │   │      (RGI)       │   │
 │   └─────────────┘   └──────────────┘   └──────────────────┘   │
 │                                                  │              │
 │                                                  ▼              │
 │                        ┌──────────────┐   ┌──────────────┐    │
 │                        │  05 · Reports│◀──│ 04 · Mobilome│    │
 │                        │  & Matrices  │   │    (MGE)     │    │
 │                        └──────────────┘   └──────────────┘    │
 │                                                                 │
 └─────────────────────────────────────────────────────────────────┘
```

---

## ⚙️ Requirements

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

## 📦 Installation

```bash
git clone https://github.com/your-repo/MetagenAMR.git
cd MetagenAMR
```

All tool environments are automatically built by Snakemake via `--use-conda` — no manual installation required.

---

## 🔧 Configuration

Edit `config/config.yaml`:

```yaml
# ── Samples ────────────────────────────────────────
samples_id:
  - HSE_2
  - HSE_4
  - HSE_5
  - HSE_6

# ── Input ──────────────────────────────────────────
reads_dir: "data/reads"          # {reads_dir}/{sample}/{sample}.fastq.gz

# ── Reference databases ────────────────────────────
kraken2_db:   "/path/to/kraken2_db"
host_genome:  "/path/to/human_genome.fa"
card_db:      "/path/to/card_db"
mobileog_db:  "/path/to/mobileOG-db_beatrix-1.6"

# ── Output settings ────────────────────────────────
fig_dpi:          300
mge_min_pident:    90
mge_min_cov:       90
```

---

## 🚀 Usage

```bash
# Complete pipeline
bash run_AMR_pipeline.sh

# Specific target
bash run_AMR_pipeline.sh -t <target>

# Dry-run (preview jobs)
bash run_AMR_pipeline.sh -n

# Custom cores
bash run_AMR_pipeline.sh -t taxonomy_all -c 16
```

### Available targets

| Target | Module | Description |
|---|---|---|
| `quality_control` | QC | Host depletion, filtering, assembly, QC metrics |
| `taxonomy_kraken` | Taxonomy | Kraken2 + Bracken classification |
| `taxonomy_analysis` | Taxonomy | Filtering, normalization, alpha diversity |
| `taxonomy_viz` | Taxonomy | HTML report |
| `taxonomy_all` | Taxonomy | Complete taxonomy module |
| `rgi_main_all` | Resistome | ARG detection (assembly) + matrices + report |
| `rgi_bwt_all` | Resistome | ARG detection (reads) + matrices + report |
| `rgi_compare_all` | Resistome | RGI main vs BWT comparison |
| `mge_all` | Mobilome | MGE detection + R analysis + report |
| `cooccurrence_all` | Co-occurrence | ARG × MGE contig-level co-localization |
| `all` | Full | Complete pipeline |

---

## 📁 Outputs

```
results/
├── {sample}/
│   ├── qc/                      # NanoPlot QC reports
│   ├── assembly/                # Flye assembled contigs
│   ├── rgi/                     # RGI main hits (.txt)
│   ├── rgi_bwt/                 # RGI BWT hits
│   └── mge/                     # mobileOG hits (.csv)
│
├── taxonomy/                    # Bracken matrices (S/G/F/O/C/P)
│   ├── abundance_matrix_S.tsv   # Species-level counts
│   ├── alpha_diversity.tsv      # Shannon, Simpson, Chao1, Pielou
│   └── ...
│
├── r_analysis/
│   ├── rgi/                     # ARG normalized matrices (11 files)
│   ├── rgi_bwt/                 # ARG BWT matrices
│   ├── mge/                     # MGE normalized matrices (5 files)
│   └── cooccurrence/            # ARG × MGE co-occurrence (4 files)
│
├── figures/
│   ├── taxonomy/taxonomy_report.html
│   ├── rgi/viz_rgi.html
│   ├── mge/mge_report.html
│   └── cooccurrence/cooccurrence_report.html
│
└── stats/
    └── sequencing_stats.tsv     # Per-sample QC summary table
```

---

## 🛠 Tools & Dependencies

### Bioinformatics tools

| Tool | Version | Usage | Reference |
|---|---|---|---|
| NanoPlot | ≥ 1.40 | QC metrics & read stats | De Coster et al. 2018 |
| Minimap2 | ≥ 2.24 | Host read depletion | Li 2018 |
| fastp-long | ≥ 0.23 | Adapter trimming & QC | Chen et al. 2018 |
| Flye | ≥ 2.9 | De novo assembly | Kolmogorov et al. 2019 |
| Kraken2 | ≥ 2.1 | Taxonomic classification | Wood et al. 2019 |
| Bracken | ≥ 2.7 | Abundance re-estimation | Lu et al. 2017 |
| RGI / CARD | ≥ 6.0 | ARG detection | Alcock et al. 2023 |
| DIAMOND | ≥ 2.1 | MGE alignment | Buchfink et al. 2021 |
| Prodigal | ≥ 2.6 | ORF prediction | Hyatt et al. 2010 |
| mobileOG-db | beatrix-1.6 | MGE database | Dong et al. 2022 |

### R & Python packages

| Package | Usage |
|---|---|
| `phyloseq` v1.44 | Alpha diversity framework |
| `vegan` v2.6-4 | Diversity indices (Shannon, Simpson, Chao1, Pielou) |
| `tidyverse` | Data wrangling |
| `pandas` / `numpy` | Matrix processing |
| `matplotlib` / `seaborn` | Figure generation |
| `scikit-learn` | CLR + PCA (beta diversity) |

---

## 📖 Citation

If you use MetagenAMR in your research, please cite:

> Gnimadi TAC et al. (2026) Manuscript in preparation.

**Please also cite key dependencies:**

```
McMurdie & Holmes (2013) phyloseq. PLoS ONE 8(4):e61217
Oksanen et al. (2022) vegan: Community Ecology Package. R package v2.6-4
Alcock et al. (2023) CARD. Nucleic Acids Research 51(D1):D690–D699
Dong et al. (2022) mobileOG-db. Applied and Environmental Microbiology
Kolmogorov et al. (2019) Flye. Nature Biotechnology 37:540–546
McMurdie & Holmes (2014) Waste not, want not. PLoS Comput Biol 10(4):e1003531
```

---

## 📄 License

```
MIT License — Copyright (c) 2026 Thibaut Armel Chérif GNIMADI et al.
```

See [LICENSE](LICENSE) for full terms.

---

<div align="center">
<sub>MetagenAMR · Snakemake · Conda · Nanopore R10.4 · CARD · mobileOG-db beatrix-1.6</sub>
</div>
