# Figure 3 - Impaired APRIL signaling in normal plasma cells

Code to reproduce Figure 3 from the integrated peripheral-blood scRNA-seq (1,109,633 clean cells; QC_removed, doublets, Platelets and CLL are excluded by these scripts; the deposit retains them with explicit labels for traceability).

| Panel | Script | Description |
|-------|--------|-------------|
| 3B (UMAP) | `01_Figure3B.py` | Full UMAP with Annotation_Level_2 cell-type labels |
| 3C (myeloid APRIL) | `02_Figure3C.py` | Myeloid TNFSF13 expression across HD/MGUS/SMM, pre vs post vaccination, age+sex adjusted |
| 3D (BM external validation) | `03_Figure3D.py` | Boiarsky et al. GSE193531: APRIL module in non-malignant bone marrow plasma cells across NBM/MGUS/SMM |
| 3E (tumor vs normal PCs) | `04_Figure3E.py` | SWIFT-seq (Lightbody et al., dbGaP phs003855): APRIL module in tumor vs non-malignant plasma cells, paired within sample, BM and PB |

## Statistical framework

Panel 3C uses an age- and sex-adjusted rank-based ANCOVA versus HD, BH-corrected across the six contrasts of the panel. Panel 3D uses a two-sided Jonckheere-Terpstra ordered-trend test across NBM -> MGUS -> SMM, plus two-sided Wilcoxon rank-sum cross-stage contrasts versus NBM with BH correction; the external dataset carries no age or sex metadata, so these contrasts are unadjusted. Panel 3E uses a two-sided paired Wilcoxon signed-rank test within each compartment, BH-corrected across the two compartments. Significance threshold q (or p) < 0.1.

## Inputs

- `scRNAseq_IMPACT_Zenodo.h5ad` — the de-identified comprehensive scRNA-seq object (1.4M cells, 42k genes, normalized log1p X + int64 `counts` layer + `obsm['X_umap']` + obs columns `Annotation_Level_1`, `Annotation_Level_2`, `Diagnosis`, `Timepoint`, `TreatmentStatus`). Download from Zenodo (DOI [10.5281/zenodo.18989222](https://doi.org/10.5281/zenodo.18989222)) and point `SCRNA_DIR` in `../config.py` at the containing directory. If you mirror the deposit to an internal Google Cloud Storage bucket, `gsutil cp gs://your-bucket/impact_data/scRNAseq_IMPACT_Zenodo.h5ad .` puts it where the scripts expect.
- External validation (panel 3D): `data/external/GSE193531_cell-level-metadata.csv` (companion cell-level malignant/normal calls from inferCNV) shipped in the Zenodo deposit; the UMI matrix is fetched from GSE193531 directly.
- Panel 3E: `data/external/swiftseq_april_persample_deid.csv` -- de-identified per-sample summary
  of the SWIFT-seq plasma-cell cohort (Lightbody et al., Nat Cancer 2025). Carries only
  de-identified sample and patient keys, compartment, grouped disease stage, per-population cell
  counts and mean module scores. Cell-level scoring was performed upstream on the source object;
  the primary data are in dbGaP phs003855.v1.p1 (controlled access). Shipped in the Zenodo
  deposit, not in git (`data/` is gitignored).

## Run

```bash
cd 03_Figure3
python 01_Figure3B.py
python 02_Figure3C.py
python 03_Figure3D.py
python 04_Figure3E.py
```
