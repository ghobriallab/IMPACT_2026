# IMPACT: Immune Dysfunction in Myeloma Precursors

Code repository for reproducing figures in:

**"Myeloma precursors are associated with suboptimal immune responses to vaccination"**

## Overview

The IMPACT study examines immune dysfunction in myeloma precursor conditions (MGUS and SMM) through multi-omic profiling of up to 731 individuals, including longitudinal serology, single-cell RNA sequencing of ~1 million immune cells, T cell receptor sequencing, and plasma proteomics.

The scripts in this repository reproduce every main and supplementary figure panel. Throughout, the significance threshold is q (or p) < 0.1, and cross-sectional disease-vs-HD contrasts are age- and sex-adjusted via rank-based ANCOVA. Paired pre/post comparisons are inherently controlled for age and sex (each individual is their own baseline). Ordered-trend statistics across the disease continuum use the Jonckheere-Terpstra test on age+sex residuals (helper in `config.R`).

## Data Access

### Zenodo (processed data)

The processed single-cell RNA-seq object and tabular data are deposited on Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18989222.svg)](https://doi.org/10.5281/zenodo.18989222)

| File | Description | Size |
|------|-------------|------|
| `scRNAseq_IMPACT_Zenodo.h5ad` | De-identified comprehensive scRNA-seq object: 1,433,497 cells, 42,090 genes, normalized log1p `X`, `counts` layer (int64), `obsm['X_umap']`, and obs columns `Annotation_Level_1`, `Annotation_Level_2`, `Diagnosis`, `Timepoint`, `TreatmentStatus`. Figure 3 / 5C scripts filter to 1,109,633 clean cells (the deposit retains doublets/QC-failed/Platelets/CLL with explicit labels so the cohort is traceable). | ~9.6 GB |
| `subclusters/` | Five de-identified per-lineage subcluster objects (B-cell, T-cell, monocyte, dendritic-cell, NK) used by Supplementary Figure 2. | ~100 MB total |
| Tabular files (`elisa/`, `olink/`, `tcr/`, `metadata/`, `external/`, plus `il1b_response_genes_human.csv`, `hvg_2678_genes.txt`, etc.) | De-identified inputs for all serology, Olink and TCR panels. | tens of MB |

After downloading, edit `SCRNA_DIR` in `config.py` and `config.R` to point to the directory containing the files, and place the tabular files under a local `data/` folder matching the layout below.

If you mirror the deposit to an internal Google Cloud Storage bucket, the following one-liners reproduce the expected layout:

```bash
# scRNA-seq h5ad + subcluster objects
gsutil cp gs://your-bucket/impact_data/scRNAseq_IMPACT_Zenodo.h5ad /path/to/scrnaseq_data/
gsutil cp -r gs://your-bucket/impact_data/subclusters /path/to/scrnaseq_data/subclusters

# Tabular inputs
gsutil cp -r gs://your-bucket/impact_data/ data/
```

### GEO (raw data)

Raw sequencing data are deposited at GEO under accession GSEXXXXXX.

### External validation datasets

| Accession | Reference | Used in |
|-----------|-----------|---------|
| [GSE193531](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193531) | Boiarsky et al., *Nat Commun* 2022 | Figure 3F (= manuscript Figure 3E) |
| [GSE124310](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124310) | Zavidij et al., *Nat Cancer* 2020 | Supplementary Figure 3 |
| [GSE205101 / GSE173644](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173644) | Stephenson et al., *J Immunol* 2022 | Supplementary Figure 4 |

Cell-level metadata for GSE193531 and the GSE173644 time-course matrix are included in the Zenodo deposit; GSE124310 sample matrices are downloaded directly by `06_SupplementalFigures/03_SupFig3.py`.

### Tabular data layout

```
data/
├── elisa/
│   ├── elisa_spike_post2nd.csv
│   ├── elisa_spike_post3rd.csv
│   ├── elisa_serial_titers_all.csv
│   ├── elisa_serial_titers_filtered.csv
│   ├── elisa_mmr.csv
│   └── elispot_spike_cef.csv
├── olink/
│   ├── olink_paired_prepost.csv
│   ├── olink_cytokines.csv
│   └── olink_summary_tumor_burden.csv
├── tcr/
│   └── tcr_clonotype_proportions.rds
├── metadata/
│   └── metadata_deidentified.csv
├── external/
│   ├── GSE193531_cell-level-metadata.csv
│   ├── PrePostTEC_olink_deid.csv
│   └── PrePostTEC_cohort_deid.csv
├── il1b_response_genes_human.csv
├── hvg_2678_genes.txt
├── smm_risk_strat.csv
└── suppfig5_risk_genetics_titer.csv
```

## Setup

### Python

```bash
conda create -n impact python=3.10
conda activate impact
pip install -r requirements.txt
```

Pinned key versions: `scanpy==1.11.5`, `anndata`, `pandas`, `numpy`, `scipy`, `statsmodels`, `matplotlib`, `seaborn`.

### R

```r
source("environment.R")
```

Installs the union of CRAN packages used across the R figure scripts.

### Configuration

Edit the path constants in `config.py` and `config.R` to point at your local data directory:

```python
# config.py
SCRNA_DIR = Path("/path/to/scrnaseq_data")
```

```r
# config.R
SCRNA_DIR <- "/path/to/scrnaseq_data"
```

## Figure Scripts

All scripts output figures to the gitignored `figures/` directory at the repo root.

Folders carry a numeric prefix matching their manuscript-figure order, and scripts inside each folder carry a numeric prefix matching their panel order. The original `FigureNX.R` / `FigureNX.py` panel name is preserved in the file name so the manuscript-to-script mapping stays explicit. A panel-letter caveat applies to Figure 3: the script names are historical, and the table below maps each script to the panel it produces in the final manuscript.

### Main Figures

| Manuscript panel | Script | Language | Description |
|------------------|--------|----------|-------------|
| 1B | `01_Figure1/01_Figure1B.R` | R | SARS-CoV-2 spike IgG titers 2-8 weeks post-2nd dose |
| 1C | `01_Figure1/02_Figure1C.R` | R | SARS-CoV-2 spike IgG titers 2-4 months post-2nd dose |
| 1D | `01_Figure1/03_Figure1D.R` | R | Per-individual serial titer trajectories |
| 1E | `01_Figure1/04_Figure1E.R` | R | Linear mixed-effects waning slope (MGUS vs SMM-untreated) |
| 1F | `01_Figure1/05_Figure1F.R` | R | Spike IgG titers after 3rd dose (SMM split into Untreated / Treated) |
| 2B | `02_Figure2/01_Figure2B.R` | R | MMR antibody titers (HD vs SMM) |
| 2D | `02_Figure2/02_Figure2D.R` | R | Plasma APRIL (TNFSF13) Olink levels, paired pre/post |
| 2E | `02_Figure2/03_Figure2E.R` | R | Post-vaccination APRIL vs serum M-spike (Spearman, treatment-naive SMM) |
| 2F | `02_Figure2/04_Figure2F.py` | Python | Paired pre/post-teclistamab plasma Olink for APRIL / BAFF / sBCMA (n=10 HRSMM on the teclistamab arm of Immuno-PRISM, NCT05469893); composite 3-protein decoy-sink readout |
| 3B (UMAP) | `03_Figure3/01_Figure3C.py` | Python | Full UMAP with Annotation_Level_2 labels |
| 3C (myeloid APRIL) | `03_Figure3/02_Figure3D.py` | Python | Myeloid TNFSF13 expression, pre vs post, age+sex adjusted |
| 3D (B-cell APRIL module) | `03_Figure3/03_Figure3E.py` | Python | APRIL-responsive module score, treatment-naive HD/MGUS/SMM |
| 3E (BM external validation) | `03_Figure3/04_Figure3F.py` | Python | Boiarsky GSE193531: APRIL module in normal vs malignant PCs across NBM/MGUS/SMM/NDMM |
| 4B | `04_Figure4/01_Figure4B.R` | R | Spike-specific TCR clonotypes (paired pre vs post) |
| 4C | `04_Figure4/02_Figure4C.R` | R | CEF-specific TCR clonotypes (recall control) |
| 4E | `04_Figure4/03_Figure4E.R` | R | IFN-gamma ELISPOT (HD vs SMM, age+sex adjusted) |
| 5A / 5B | `05_Figure5/01_Figure5AB.R` | R | Paired IL-1B and IL-18 Olink levels |
| 5C | `05_Figure5/02_Figure5C.py` | Python | IL-1B response gene signature score (scRNA-seq) |
| 5D / 5E / 5F | `05_Figure5/03_Figure5DEF.R` | R | Paired DDX58, NUB1 and MMP7 Olink levels |

### Supplementary Figures

| Manuscript panel | Script | Language | Description |
|------------------|--------|----------|-------------|
| S1 | `06_SupplementalFigures/01_SupFig1.py` | Python | Vaccine response by SMM 20/2/20 risk tier (LR / IR / HR), with per-tier treated/treatment-naive n annotated |
| S2 | `06_SupplementalFigures/02_SupFig2.py` | Python | Per-lineage cell-type annotation UMAPs + canonical marker-gene heatmaps |
| S3 | `06_SupplementalFigures/03_SupFig3.py` | Python | Bone-marrow myeloid TNFSF13 (APRIL) across HD/MGUS/SMM/MM (Zavidij GSE124310) |
| S4 | `06_SupplementalFigures/04_SupFig4.py` | Python | External validation of the APRIL-responsive gene signature (GSE205101 / GSE173644) |
| S5 | `06_SupplementalFigures/05_SupFig5.py` | Python | Individual APRIL-responsive gene violins (HD vs SMM, post-vaccination) |
| S6 | `06_SupplementalFigures/06_SupFig6.R` | R | Sample shipping (FedEx vs not-shipped) does not confound SARS-CoV-2 spike-specific TCR clonotype frequencies (paired-design control for Figure 4B) |

### Execution order

Most scripts are independent. The only intra-folder dependency is:

- **`03_Figure5DEF.R`** sources **`01_Figure5AB.R`** for its data prep and the `plot_cytokine()` helper.

Scripts requiring the scRNA-seq h5ad file: Figure 3 (all four scripts), Figure 5C, and Supplementary Figures 2 and 5.

## Reproducibility

This repository carries an automated reproducibility-check workflow at `.github/workflows/run-reproducibility-check.yml`. It runs `tools/check_reproducibility.py` weekly, scores the repo on six categories (step ordering, documentation, path hygiene, GCS data handling, naming conventions, PHI/credential safety), and appends the result to `reproducibility.md`. Locally:

```bash
python3 tools/check_reproducibility.py .
```

## Citation

If you use this code or data, please cite:
