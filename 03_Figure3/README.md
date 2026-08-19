# Figure 3: Impaired APRIL signaling in plasma cells

Panels 3B and 3C come from the integrated peripheral-blood scRNA-seq (1,109,633 cells after
excluding QC-failed cells, doublets, platelets and CLL; the deposited object retains these under
explicit labels so the cohort stays traceable). Panels 3D and 3E come from the SWIFT-seq
bone-marrow and peripheral-blood plasma-cell cohort (Lightbody et al., Nat Cancer 2025).

| Panel | Script | Description |
|-------|--------|-------------|
| 3B | `01_Figure3B.py` | Full UMAP with Annotation_Level_2 cell-type labels |
| 3C | `02_Figure3C.py` | Myeloid TNFSF13 expression across HD/MGUS/SMM, pre vs post vaccination |
| 3D | `03_Figure3D.py` | APRIL-responsive module in non-malignant plasma cells across NBM/MGUS/SMM/NDMM, bone marrow and peripheral blood |
| 3E | `04_Figure3E.py` | APRIL-responsive module in tumor vs non-malignant plasma cells, paired within sample |

Panels 3D and 3E keep one sample per participant per compartment, the earliest serial timepoint,
so that repeatedly sampled participants are not counted more than once.

## Statistics

Panel 3C uses age- and sex-adjusted rank-based ANCOVA against HD, Benjamini-Hochberg corrected
across the six contrasts of the panel. Panel 3D uses a two-sided Jonckheere-Terpstra ordered-trend
test across NBM to NDMM, plus two-sided Wilcoxon rank-sum contrasts against NBM with BH correction
within each compartment; the external cohort carries no age or sex metadata, so these are
unadjusted. Panel 3E uses a two-sided paired Wilcoxon signed-rank test within each compartment, BH
corrected across the two compartments. Significance threshold q (or p) < 0.1.

## Inputs

- `scRNAseq_IMPACT_Zenodo.h5ad` for panels 3B and 3C. Download it from Zenodo
  (DOI [10.5281/zenodo.18989222](https://doi.org/10.5281/zenodo.18989222)) and point `SCRNA_DIR`
  in `../config.py` at the containing directory.
- `data/external/swiftseq_april_persample_deid.csv` for panels 3D and 3E: a de-identified
  per-sample summary carrying sample and patient keys, compartment, disease stage, serial
  timepoint index, per-population cell counts and mean module scores. Cell-level scoring was done
  upstream on the source object with `sc.tl.score_genes()` and the gene list in Supplementary
  Table 5; the primary data are in dbGaP phs003855.v1.p1 (controlled access). The summary ships in
  the Zenodo deposit.

## Run

```bash
cd 03_Figure3
python3 01_Figure3B.py    # writes ../figures/Figure3B.{png,pdf,svg}
python3 02_Figure3C.py
python3 03_Figure3D.py
python3 04_Figure3E.py
```
