# Supplementary Figures

Numbered in the order they are first cited in the manuscript.

| Panel | Script | Description |
|-------|--------|-------------|
| S1 | `01_SupFig1.py` | Vaccine response by SMM 20/2/20 risk tier (LR / IR / HR), with the treated and treatment-naive n annotated per tier |
| S2 | `02_SupFig2.py` | Per-lineage cell-type annotation UMAPs and canonical marker-gene heatmaps |
| S3A | `03_SupFig3A.py` | Myeloid TNFSF13 (APRIL) by shipment status in SMM, the shipment effect itself |
| S3B | `03_SupFig3B.py` | Figure 3C reproduced on shipped samples only, holding shipment constant |
| S3C | `03_SupFig3C.py` | Bone-marrow myeloid TNFSF13 (APRIL) across HD/MGUS/SMM/MM (Zavidij et al. 2020, GSE124310) |
| S4A | `04_SupFig4A.py` | APRIL-responsive signature in APRIL-stimulated plasmablasts (GSE173644) |
| S4B | `04_SupFig4B.py` | APRIL-responsive module in non-malignant bone-marrow plasma cells across NBM/MGUS/SMM (Boiarsky et al. 2022, GSE193531) |
| S5 | `05_SupFig5.R` | Shipment control for Figure 4B: spike-specific TCR clonotype frequencies by shipment status |
| S6 | `06_SupFig6.py` | Shipment control for Figure 5C: the IL-1B response signature by shipment status |

## Statistics

The same conventions as the main figures: age- and sex-adjusted rank-based ANCOVA for
cross-sectional disease-vs-HD contrasts where age and sex are available, Jonckheere-Terpstra for
ordered-trend claims, Wilcoxon signed-rank for paired comparisons, and q (or p) < 0.1 throughout.
Where an external dataset carries no age or sex metadata (S3C, S4A, S4B), contrasts are unadjusted
two-sided Wilcoxon rank-sum tests and each script docstring says so.

## Shipment controls

S3A, S3B, S5 and S6 all read shipment status from `Shipment_FedEx` in Supplementary Table 4, where 1 is
shipped and 0 is not. A blank field is excluded rather than assumed shipped, so a missing record
never counts as evidence either way.

S3A shows the shipment effect directly, shipped versus not shipped within SMM, which is the only
Figure 3C group where shipment varies: every healthy donor and all but one participant with MGUS
was shipped. S3B then reproduces Figure 3C using only samples documented as shipped, holding
shipment constant, and the disease comparison stays null.
S6 repeats the Figure 5C paired test within each shipment stratum, then asks separately whether
the pre-to-post change depends on shipment; both panels are restricted to treatment-naive SMM,
which is the only Figure 5C group containing both strata, since every paired HD and MGUS
participant was shipped.

S3B and S6 each assert that they reproduce their parent panel before reporting any stratified
result, and S3A shares S3B's aggregation code so the two panels cannot describe different data.
S5 reports an unadjusted two-sided Wilcoxon rank-sum test at each timepoint, with clonotype
proportions as percentages to match the Figure 4B axis.

## Inputs

| File | Panels |
|------|--------|
| `data/elisa/elisa_spike_post2nd.csv`, `elisa_spike_post3rd.csv`, `data/smm_risk_strat.csv` | S1 |
| `scRNAseq_IMPACT_Zenodo.h5ad` and the per-lineage objects in `SUBCLUSTER_DIR` | S2 |
| `data/metadata/Supplementary_Table_4_scRNAseq_sample_list.csv` and `scRNAseq_IMPACT_Zenodo.h5ad` | S3A, S3B, S6 |
| `data/external/zavidij_bm/matrices/` (GSE124310, downloaded from GEO) | S3C |
| `data/external/GSE173644_timecourse.txt.gz` | S4A |
| `data/external/GSE193531_umi-count-matrix.csv.gz` and `GSE193531_cell-level-metadata.csv` | S4B |
| `data/tcr/tcr_clonotype_proportions.rds`, `data/tcr/tcr_shipping_status.csv` | S5 |
| `data/il1b_response_genes_human.csv`, `data/hvg_2678_genes.txt` | S6 |

S6 imports the cell filter, effect size and statistics helpers from `05_Figure5/02_Figure5C.py`,
so the two figures cannot drift apart. Every input except the GSE124310 matrices ships in the
Zenodo deposit (DOI [10.5281/zenodo.18989222](https://doi.org/10.5281/zenodo.18989222)).

## Run

```bash
cd 06_SupplementalFigures
python3 01_SupFig1.py     # writes ../figures/SupFig1.{png,pdf,svg}
python3 02_SupFig2.py
python3 03_SupFig3A.py
python3 03_SupFig3B.py
python3 03_SupFig3C.py
python3 04_SupFig4A.py
python3 04_SupFig4B.py
Rscript 05_SupFig5.R
python3 06_SupFig6.py
```
