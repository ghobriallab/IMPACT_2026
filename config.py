"""IMPACT — Configuration. Edit these paths to match your local environment.

Purpose:      Shared configuration for Python figure scripts: paths to the de-identified scRNA-seq h5ad files, data and figures directories, and the manuscript color palette.

Inputs:       Edit SCRNA_DIR to point to the folder containing scRNAseq_IMPACT_Zenodo.h5ad (downloaded from Zenodo DOI 10.5281/zenodo.18989223).

Outputs:      Module-level constants imported by every Python figure script (SCRNA_DIR, H5AD_*, DATA_DIR, FIGURES_DIR, COLORS).

Dependencies: Python 3.10+; pathlib (stdlib).
"""
from pathlib import Path

# Base directory for this repository (auto-detected from this file's location)
REPO_DIR = Path(__file__).parent

# Base directory containing scRNA-seq h5ad files
# Download scRNAseq_IMPACT_Zenodo.h5ad from Zenodo (https://doi.org/10.5281/zenodo.18989223)
# and edit this path to the directory containing it.
SCRNA_DIR = Path("/path/to/scrnaseq_data")

# Specific h5ad files
# consolidated to ONE comprehensive de-identified object. scRNAseq_IMPACT_Zenodo.h5ad carries
# all 1,433,497 cells (doublets/QC retained but LABELED in Annotation_Level_2; CLL removed), 42,090 genes,
# normalized log1p X + int64 'counts' layer + obsm['X_umap'] + obs (Annotation_Level_1/2, Diagnosis,
# Timepoint, TreatmentStatus). Fig3 scripts filter to clean cells (-> 1,109,633); Figure5C re-normalizes
# from counts and subsets to data/hvg_2678_genes.txt before score_genes. All handles point at this file.
H5AD_CELLXGENE = SCRNA_DIR / "scRNAseq_IMPACT_Zenodo.h5ad"
H5AD_NORM = SCRNA_DIR / "scRNAseq_IMPACT_Zenodo.h5ad"
H5AD_ANNOTATED = SCRNA_DIR / "scRNAseq_IMPACT_Zenodo.h5ad"
H5AD_IL1B = SCRNA_DIR / "scRNAseq_IMPACT_Zenodo.h5ad"  # Figure5C subsets to the 2,678-HVG list (control pool)
SUBCLUSTER_DIR = SCRNA_DIR / "subclusters"  # SupFig1 only; per-lineage subcluster objects (see README/notes)

# Data files (download from Zenodo: https://doi.org/10.5281/zenodo.18989223)
DATA_DIR = REPO_DIR / "data"

# Output directory for figures
FIGURES_DIR = REPO_DIR / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

# Metadata
METADATA_CSV = DATA_DIR / "metadata" / "metadata_deidentified.csv"

# Color palette (consistent across all figures)
COLORS = {
    "HD": "#4DBBD5",
    "MGUS": "#F39B7F",
    "SMM": "#E64B35",
    "SMM (Untreated)": "#E64B35",
    "SMM (Treated)": "#91D1C2",
    "MM": "#3C5488",
}
