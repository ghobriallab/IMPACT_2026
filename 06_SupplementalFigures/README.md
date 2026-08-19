# Supplementary Figures

Code to reproduce the Supplementary Figures. Supplementary figures are numbered in citation order in the manuscript.

| Supp. Fig. | Script | Description |
|------------|--------|-------------|
| S1 | `01_SupFig1.py` | Vaccine response by SMM 20/2/20 risk tier (Fig 1B/1C/1F equivalents stratified into LR / IR / HR-SMM), with per-tier treated-vs-naive n's annotated on each column |
| S2 | `02_SupFig2.py` | Per-lineage cell-type annotation UMAPs + canonical marker-gene heatmaps for the integrated scRNA-seq object |
| S3A | `03_SupFig3A_shipment.py` | Figure 3C reproduced on shipped samples only; shipment control for Figure 3C |
| S3B | `03_SupFig3.py` | Bone-marrow myeloid TNFSF13 (APRIL) expression across HD/MGUS/SMM/MM (Zavidij et al. 2020, GSE124310) |
| S4A | `04_SupFig4A.py` | External validation of the APRIL-responsive gene signature in APRIL-stimulated plasmablasts (GSE173644) |
| S4B | `04_SupFig4B.py` | APRIL-responsive module in non-malignant bone-marrow plasma cells across NBM/MGUS/SMM (Boiarsky GSE193531); was main Figure 3D |
| S5 | `05_SupFig5.R` | Sample shipping (shipped vs not shipped) does not confound SARS-CoV-2 spike-specific TCR clonotype frequencies; paired-design control for Figure 4B |
| S6 | `06_SupFig6.py` | Sample shipment and the IL-1beta response signature; shipment control for Figure 5C |

## Statistical framework

Same conventions as the main figures: age- and sex-adjusted rank-based ANCOVA for cross-sectional disease-vs-HD contrasts (where age and sex metadata are available); Jonckheere-Terpstra ordered-trend test for monotonic-gradient claims; Wilcoxon signed-rank for paired comparisons. Significance threshold q (or p) < 0.1. Where an external dataset lacks age and sex metadata (S3, Zavidij), the comparison is an unadjusted two-sided Wilcoxon rank-sum test and the script docstring says so.

S5 reports an unadjusted two-sided Wilcoxon rank-sum test at each timepoint. Clonotype proportions are expressed as percentages, matching the y-axis of Figure 4B.

The three shipment controls (S3A, S5, S6) all define shipment from `Shipment_FedEx` in Supplementary Table 4, where 1 is shipped and 0 is not shipped. A blank field is excluded rather than assumed shipped, so a missing record never counts as evidence either way. S3A reproduces Figure 3C using only samples documented as shipped, holding shipment constant; S6 repeats the Figure 5C paired test within each shipment stratum and asks separately whether the pre-to-post change depends on shipment, both restricted to treatment-naive SMM because every paired HD and MGUS participant was shipped. Both assert that they reproduce the published panel before reporting any stratified result.

## Inputs

- `scRNAseq_IMPACT_Zenodo.h5ad` and the per-lineage subcluster objects in `SUBCLUSTER_DIR`, for the lineage UMAPs and marker heatmaps (S2). S2 is the only supplementary figure that reads the single-cell object.
- `data/elisa/elisa_spike_post2nd.csv`, `data/elisa/elisa_spike_post3rd.csv` and `data/smm_risk_strat.csv` for the risk-tier panels (S1).
- External GSE124310 (Zavidij) sample matrices for S3B; GSE173644 in-vitro APRIL stimulation time course for S4A and GSE193531 (Boiarsky) cell-level metadata and count matrix for S4B, both shipped in the Zenodo deposit.
- `data/tcr/tcr_clonotype_proportions.rds` and `data/tcr/tcr_shipping_status.csv` for S5.
- `data/metadata/Supplementary_Table_4_scRNAseq_sample_list.csv` for the `Shipment_FedEx` flag and the age/sex covariates used by S3A and S6, and `scRNAseq_IMPACT_Zenodo.h5ad` for their expression values. S6 additionally reads `data/il1b_response_genes_human.csv` and `data/hvg_2678_genes.txt`, and imports the cell filter, effect size and statistics helpers from `05_Figure5/02_Figure5C.py` so the two figures cannot drift apart.

All in-house inputs are hosted on Zenodo (DOI [10.5281/zenodo.18989222](https://doi.org/10.5281/zenodo.18989222)). To mirror the deposit to an internal Google Cloud Storage bucket, copy the whole Zenodo file tree into your bucket; `gsutil cp -r gs://your-bucket/impact_data/ data/` then reproduces the expected layout.

## Run

```bash
cd 06_SupplementalFigures
python3 01_SupFig1.py
python3 02_SupFig2.py
python3 03_SupFig3A_shipment.py
python3 03_SupFig3.py
python3 04_SupFig4A.py
python3 04_SupFig4B.py
Rscript 05_SupFig5.R
python3 06_SupFig6.py
```
