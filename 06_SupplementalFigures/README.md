# Supplementary Figures

Code to reproduce the Supplementary Figures. Supplementary figures are numbered in citation order in the manuscript.

| Supp. Fig. | Script | Description |
|------------|--------|-------------|
| S1 | `SupFig1.py` | Vaccine response by SMM 20/2/20 risk tier (Fig 1B/1C/1F equivalents stratified into LR / IR / HR-SMM), with per-tier treated-vs-naive n's annotated on each column |
| S2 | `SupFig2.py` | Per-lineage cell-type annotation UMAPs + canonical marker-gene heatmaps for the integrated scRNA-seq object |
| S3 | `SupFig3.py` | Bone-marrow myeloid TNFSF13 (APRIL) expression across HD/MGUS/SMM/MM (Zavidij et al. 2020, GSE124310) |
| S4 | `SupFig4.py` | External validation of the APRIL-responsive gene signature in independent B-cell scRNA-seq cohorts |
| S5 | `SupFig5.py` | Individual APRIL-responsive gene violin plots (HD vs SMM, post-vaccination) |

## Statistical framework

Same conventions as the main figures: age- and sex-adjusted rank-based ANCOVA for cross-sectional disease-vs-HD contrasts (where age and sex metadata are available); Jonckheere–Terpstra ordered-trend test for monotonic-gradient claims; Wilcoxon signed-rank for paired comparisons. Significance threshold q (or p) < 0.1. Where external datasets lack age/sex metadata (e.g. SupFig3 Zavidij), this is explicitly noted in the script docstring.

## Inputs

- `scRNAseq_IMPACT_Zenodo.h5ad` for the lineage UMAPs (S2) and the APRIL gene violins (S5).
- Per-lineage subcluster objects in `SUBCLUSTER_DIR` (S2).
- `data/elisa/elisa_spike_post2nd.csv`, `elisa_spike_post3rd.csv`, `data/smm_risk_strat.csv` for the risk-tier panels (S1).
- External GSE124310 (Zavidij) sample matrices for S3; GSE205101 / GSE173644 in-vitro APRIL stimulation matrices for S4.

All in-house inputs are hosted on Zenodo (DOI [10.5281/zenodo.18989222](https://doi.org/10.5281/zenodo.18989222)). To mirror the deposit to an internal Google Cloud Storage bucket, copy the whole Zenodo file tree into your bucket and `gsutil cp -r gs://your-bucket/impact_data/ data/` reproduces the expected layout for the SupplementalFigures scripts.

## Run

```bash
cd SupplementalFigures
python SupFig1.py
python SupFig2.py
python SupFig3.py
python SupFig4.py
python SupFig5.py
```
