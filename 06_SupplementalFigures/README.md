# Supplementary Figures

Code to reproduce the Supplementary Figures. Supplementary figures are numbered in citation order in the manuscript.

| Supp. Fig. | Script | Description |
|------------|--------|-------------|
| S1 | `01_SupFig1.py` | Vaccine response by SMM 20/2/20 risk tier (Fig 1B/1C/1F equivalents stratified into LR / IR / HR-SMM), with per-tier treated-vs-naive n's annotated on each column |
| S2 | `02_SupFig2.py` | Per-lineage cell-type annotation UMAPs + canonical marker-gene heatmaps for the integrated scRNA-seq object |
| S3 | `03_SupFig3.py` | Bone-marrow myeloid TNFSF13 (APRIL) expression across HD/MGUS/SMM/MM (Zavidij et al. 2020, GSE124310) |
| S4 | `04_SupFig4.py` | External validation of the APRIL-responsive gene signature in APRIL-stimulated plasmablasts (GSE173644) |
| S5 | `05_SupFig5.R` | Sample shipping (shipped vs not shipped) does not confound SARS-CoV-2 spike-specific TCR clonotype frequencies; paired-design control for Figure 4B |

## Statistical framework

Same conventions as the main figures: age- and sex-adjusted rank-based ANCOVA for cross-sectional disease-vs-HD contrasts (where age and sex metadata are available); Jonckheere-Terpstra ordered-trend test for monotonic-gradient claims; Wilcoxon signed-rank for paired comparisons. Significance threshold q (or p) < 0.1. Where an external dataset lacks age and sex metadata (S3, Zavidij), the comparison is an unadjusted two-sided Wilcoxon rank-sum test and the script docstring says so.

S5 reports an unadjusted two-sided Wilcoxon rank-sum test at each timepoint. Clonotype proportions are expressed as percentages, matching the y-axis of Figure 4B.

## Inputs

- `scRNAseq_IMPACT_Zenodo.h5ad` and the per-lineage subcluster objects in `SUBCLUSTER_DIR`, for the lineage UMAPs and marker heatmaps (S2). S2 is the only supplementary figure that reads the single-cell object.
- `data/elisa/elisa_spike_post2nd.csv`, `data/elisa/elisa_spike_post3rd.csv` and `data/smm_risk_strat.csv` for the risk-tier panels (S1).
- External GSE124310 (Zavidij) sample matrices for S3, downloaded by the script; GSE173644 in-vitro APRIL stimulation time course for S4, shipped in the Zenodo deposit.
- `data/tcr/tcr_clonotype_proportions.rds` and `data/tcr/tcr_shipping_status.csv` for S5.

All in-house inputs are hosted on Zenodo (DOI [10.5281/zenodo.18989222](https://doi.org/10.5281/zenodo.18989222)). To mirror the deposit to an internal Google Cloud Storage bucket, copy the whole Zenodo file tree into your bucket; `gsutil cp -r gs://your-bucket/impact_data/ data/` then reproduces the expected layout.

## Run

```bash
cd 06_SupplementalFigures
python3 01_SupFig1.py
python3 02_SupFig2.py
python3 03_SupFig3.py
python3 04_SupFig4.py
Rscript 05_SupFig5.R
```
