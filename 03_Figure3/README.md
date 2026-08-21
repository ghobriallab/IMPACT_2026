# Figure 3: Impaired APRIL signaling in plasma cells

Panels 3B and 3C come from the integrated peripheral-blood scRNA-seq (1,109,633 cells after
excluding QC-failed cells and doublets; the deposited object retains these under explicit
labels so the cohort stays traceable). Platelets and plasma cells are kept in this embedding and shown for
completeness; they are excluded from the downstream analyses. Panels 3D and 3E come from the
SWIFT-seq bone-marrow and peripheral-blood plasma-cell cohort (Lightbody et al., Nat Cancer 2025).

| Panel | Script | Description |
|-------|--------|-------------|
| 3B | `01_Figure3B.py` | Full UMAP with Annotation_Level_2 cell-type labels |
| 3C | `02_Figure3C.py` | Myeloid TNFSF13 expression across HD/MGUS/SMM, pre vs post vaccination |
| 3D | `03_Figure3D.py` | APRIL-responsive module in non-malignant plasma cells across NBM/MGUS/SMM/NDMM, bone marrow and peripheral blood |
| 3E | `04_Figure3E.py` | APRIL-responsive module in tumor vs non-malignant plasma cells, paired within sample |

Panels 3D and 3E follow the cohort's own convention for the SWIFT-seq data: baseline specimens
only, plasma cells labelled Normal, Tumor or Tumor1 (secondary tumour clones and unannotated cells
excluded), and disease stage taken from the cohort's FinalDx. A specimen is one participant and
one tissue, so the several specimens sequenced as more than one library are pooled rather than
counted twice, and no participant enters a comparison more than once per compartment. Age and sex
are carried through from the cohort metadata so that panel 3D can be adjusted for them, as every
other cross-sectional comparison in the paper is.

## Statistics

Panel 3C uses age- and sex-adjusted rank-based ANCOVA against HD, Benjamini-Hochberg corrected
across the six contrasts of the panel. Panel 3D uses a two-sided Jonckheere-Terpstra ordered-trend
test on age- and sex-adjusted residuals across NBM to NDMM, plus age- and sex-adjusted rank-based
ANCOVA contrasts against NBM with BH correction across the three contrasts of each compartment.
Panel 3E uses a two-sided paired Wilcoxon signed-rank test within each compartment, BH
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
