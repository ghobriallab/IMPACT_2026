# IMPACT: Immune Dysfunction in Myeloma Precursors

Code repository for reproducing figures in:

**"Myeloma Precursors Erode Durable Immunity"**

## Overview

The IMPACT study examines immune dysfunction in myeloma precursor conditions (MGUS and SMM) through multi-omic profiling of up to 731 individuals, including longitudinal serology, single-cell RNA sequencing of ~1 million immune cells, T cell receptor sequencing, and plasma proteomics.

## Data Access

### Zenodo (processed data)

The processed single-cell RNA-seq object and TCR data are deposited on Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18989223.svg)](https://doi.org/10.5281/zenodo.18989223)

| File | Description | Size |
|------|-------------|------|
| `scRNAseq_IMPACT_Zenodo.h5ad` | Integrated scRNA-seq (1.1M cells, 42k genes) | 18.7 GB |
| `TCR.csv` | De-identified TCR diversity and clonotype data (104 patients) | 31 KB |

After downloading `scRNAseq_IMPACT_Zenodo.h5ad`, edit the paths in `config.py` and `config.R` to point to the directory containing the file.

### GEO (raw data)

Raw sequencing data are deposited at GEO under accession GSEXXXXXX.

### External validation datasets

| Accession | Reference | Used in |
|-----------|-----------|---------|
| [GSE193531](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193531) | Boiarsky et al., *Nat Commun* 2022 | Figure 3F |
| [GSE173644](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173644) | Stephenson et al., *J Immunol* 2022 | Sup. Figure 2 |

Cell-level metadata for GSE193531 and the GSE173644 time-course matrix are included in the Zenodo deposit.

### Tabular data files (Zenodo)

All tabular data files are hosted on Zenodo alongside the scRNA-seq object. After downloading, place them in the `data/` directory to match the expected layout:

```
data/
├── elisa/
│   ├── elisa_spike_post2nd.csv
│   ├── elisa_spike_post3rd.csv
│   ├── elisa_spike_post2nd_raw.csv
│   ├── elisa_serial_titers_all.csv
│   ├── elisa_serial_titers_filtered.csv
│   ├── elisa_mmr.csv
│   └── elispot_spike_cef.csv
├── olink/
│   ├── olink_paired_prepost.csv
│   ├── olink_cytokines.csv
│   └── olink_summary_tumor_burden.csv
├── tcr/
│   ├── tcr_clonotype_proportions.rds
│   └── tcr_diversity_downsampled.rds
├── metadata/
│   └── metadata_deidentified.csv
├── external/
│   ├── GSE193531_cell-level-metadata.csv
│   └── GSE173644_timecourse.txt.gz
└── il1b_response_genes_human.csv
```

## Setup

### Python

```bash
conda create -n impact python=3.10
conda activate impact
pip install -r requirements.txt
```

### R

```r
source("environment.R")
```

### Configuration

Edit `config.py` and `config.R` to set the path to your downloaded scRNA-seq data:

```python
# config.py
SCRNA_DIR = Path("/your/path/to/scrnaseq_data")
```

```r
# config.R
SCRNA_DIR <- "/your/path/to/scrnaseq_data"
```

## Figure Scripts

All scripts output figures to the `figures/` directory.

### Main Figures

| Figure | Script | Language | Description |
|--------|--------|----------|-------------|
| 1B | `Figure1/Figure1B.R` | R | SARS-CoV-2 IgG titers 2–8 weeks post-2nd dose |
| 1C | `Figure1/Figure1C.R` | R | SARS-CoV-2 IgG titers 2–4 months post-2nd dose |
| 1D | `Figure1/Figure1D.R` | R | Serial antibody titer trajectories |
| 1E | `Figure1/Figure1E.R` | R | Antibody waning rate (mixed-effects model) |
| 1F | `Figure1/Figure1F.R` | R | SARS-CoV-2 IgG titers after 3rd dose |
| 2B | `Figure2/Figure2B.R` | R | MMR antibody titers |
| 2D | `Figure2/Figure2D.R` | R | APRIL protein levels (Olink) |
| 2E | `Figure2/Figure2E.R` | R | APRIL vs M-spike correlation |
| 3C | `Figure3/Figure3C.py` | Python | Full UMAP embedding with cell type annotations |
| 3D | `Figure3/Figure3D.py` | Python | APRIL expression in myeloid cells |
| 3E | `Figure3/Figure3E.py` | Python | APRIL-responsive gene module scores |
| 3F | `Figure3/Figure3F.py` | Python | External validation of APRIL targets (GSE193531) |
| 4B | `Figure4/Figure4B.R` | R | Spike-specific TCR clonotypes (pre vs post) |
| 4C | `Figure4/Figure4C.R` | R | CEF-specific TCR clonotypes (pre vs post) |
| 4E | `Figure4/Figure4E.R` | R | IFN-γ ELISPOT |
| 5A–B | `Figure5/Figure5AB.R` | R | IL-1β and IL-18 protein levels (Olink) |
| 5C | `Figure5/Figure5C.py` | Python | IL-1β response gene signature (scRNA-seq) |
| 5D–F | `Figure5/Figure5DEF.R` | R | DDX58, NUB1, MMP7 protein levels (Olink) |

### Supplemental Figures

| Figure | Script | Language | Description |
|--------|--------|----------|-------------|
| S1 | `SupplementalFigures/SupFig1.py` | Python | Cell lineage UMAPs and marker heatmaps |
| S2 | `SupplementalFigures/SupFig2.py` | Python | APRIL signature validation time course (GSE173644) |
| S3 | `SupplementalFigures/SupFig3.py` | Python | Individual APRIL gene expression (HD vs SMM) |

### Execution order

Most scripts are independent. The only dependency is:

- **Figure5DEF.R** sources **Figure5AB.R** for shared data loading and the `plot_cytokine()` function.

Scripts requiring the scRNA-seq h5ad file: Figure 3C, 3D, 3E, 3F, 5C, SupFig 1, 2, 3.

## Citation

If you use this code or data, please cite:
