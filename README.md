# SenPred_HDF

[![R version](https://img.shields.io/badge/R-4.2.2-276DC3?style=flat&logo=r)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-4.3.0.1-blue?style=flat)](https://satijalab.org/seurat/)
[![scPred](https://img.shields.io/badge/scPred-1.9.2-green?style=flat)](https://github.com/powellgenomicslab/scPred)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/Published-Genome%20Medicine-red?style=flat)](https://genomemedicine.biomedcentral.com)

**Machine learning classifier to detect senescent cells in vivo from single-cell RNA-seq data.**

Built from scRNA-seq of replicative senescent and early proliferative human dermal fibroblasts (HDFs). Models are validated across 5 independent external datasets spanning skin and lung tissue.

> 📄 **Published in *Genome Medicine*** · Blizard Institute Paper of the Year 2025  
> **Authors:** Beth Hughes, Deborah Milligan, Cleo Bishop — Queen Mary University of London

---

## Overview

Cellular senescence plays a key role in ageing and disease, but identifying senescent cells *in vivo* is challenging due to the lack of single definitive markers. SenPred addresses this by training a multi-class ML classifier (Early Proliferative / Early Senescent / Deep Senescent) on well-characterised scRNA-seq data, then applying it to external datasets.

```
scRNA-seq (HDFs)
      │
      ▼
┌─────────────────────────────┐
│  QC · Filtering · Normalise  │  Scripts 1–2
└──────────────┬──────────────┘
               │
               ▼
┌─────────────────────────────┐
│   Clustering · Annotation    │  Script 3
└──────────────┬──────────────┘
               │
               ▼
┌─────────────────────────────┐
│  scPred ML Model Training    │  Scripts 4, 6, 11
│  (EP / ES / DS classes)      │
└──────────────┬──────────────┘
               │
      ┌────────┴────────┐
      ▼                 ▼
 2D Models           3D Models
 Scripts 5–10        Scripts 12–14
      │                 │
      └────────┬────────┘
               ▼
┌─────────────────────────────┐
│  Validation on external      │
│  in vivo skin + lung data    │  Scripts 15–22
│  (Tabib, Solé-Boldo,         │
│   Ganier, HCA Lung Atlas)    │
└─────────────────────────────┘
```

---

## Input data

Input data is available on GEO: **[GSE282425](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282425)**

Download and place contents into a folder called `Inputfiles/` within the R project before running.

---

## Installation

Package versions are stored in `renv.lock`. To restore the environment:

```r
# Install renv if needed
install.packages("renv")

# Restore all packages at correct versions
renv::restore()
```

---

## Usage

Run scripts in numbered order (1–22). Helper scripts are provided:

```r
# Clean/create the output folder first
source("R_scripts/95_Make_Clean.rmd")

# Then run the full pipeline
source("R_scripts/99_Run_All.rmd")
```

> ⚠️ The train/test split is random by default. For reproducible outputs, save the Seurat objects from `Output/Seurat_Objects/` and uncomment the `ReadRDS` lines in scripts 4, 6, and 9.

---

## Pipeline scripts

| Script | Description |
|--------|-------------|
| `1_HDF_data.Rmd` | Load data → Seurat object |
| `2_HDF_filtering_normalisation.Rmd` | QC, filtering, normalisation |
| `3_HDF_2DCells.Rmd` | 2D cell selection, clustering, DE |
| `4_Scpred_2D.Rmd` | Build 2D ML model (EP/DS) |
| `5_ScPred2d_Chan.Rmd` | Validate on Chan et al. dataset |
| `6_ScPred2D_Chan_EPESDS.Rmd` | Extend model with Early Senescent class |
| `7_Teo.Rmd` | Test on OIS/paracrine senescence (Teo et al.) |
| `8_Tabib.Rmd` | Validate on whole-skin dataset (Tabib et al.) |
| `9_soleboldo.Rmd` | Validate on whole-skin dataset (Solé-Boldo et al.) |
| `10_soleboldo_tabib_combined2D.Rmd` | Combine donor predictions — 2D model |
| `11_ScPred_3D.Rmd` | Build 3D ML model |
| `12–14` | Validate 3D model on Tabib, Solé-Boldo datasets |
| `15_matrisome.Rmd` | ECM/matrisome expression comparison |
| `16–18_Ganier.Rmd` | Validate + age-trend analysis (Ganier et al.) |
| `19–20` | Test on all skin cell types |
| `21–22_lungs.Rmd` | Apply SenPred to lung fibroblasts (HCA Lung Atlas) |

---

## External datasets used for validation

| Dataset | Tissue | Reference |
|---------|--------|-----------|
| Chan et al. 2022 | HDFs (multiomics) | *eLife* |
| Teo et al. 2019 | HDFs (OIS, paracrine) | *Cell Reports* |
| Tabib et al. 2018 | Whole skin | *J Invest Dermatol* |
| Solé-Boldo et al. 2020 | Whole skin (ageing) | *Commun Biology* |
| Ganier et al. 2024 | Skin + BCC | *PNAS* |
| Sikkema et al. 2023 | Lung (HCA) | *Nature Medicine* |

---

## Funding

BBSRC iCASE Unilever PhD studentship (BB/V509668/1)

---

## Contact

**Beth Hughes** · b.k.hughes@qmul.ac.uk · Queen Mary University of London  
**PI: Prof. Cleo Bishop** · c.l.bishop@qmul.ac.uk
