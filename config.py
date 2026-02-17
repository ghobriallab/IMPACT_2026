"""IMPACT 2026 — Configuration. Edit these paths to match your local environment."""

from pathlib import Path

# Base directory for this repository (auto-detected from this file's location)
REPO_DIR = Path(__file__).parent

# Base directory containing scRNA-seq h5ad files (from GEO/Zenodo)
SCRNA_DIR = Path("/path/to/scrnaseq_data")

# Specific h5ad files
H5AD_CELLXGENE = SCRNA_DIR / "scRNAseq_IMPACT_Zenodo.h5ad"
H5AD_NORM = SCRNA_DIR / "ad_all_norm_log1p.h5ad"
H5AD_ANNOTATED = SCRNA_DIR / "ad_all_hmy_pt_umap_annotated.h5ad"
SUBCLUSTER_DIR = SCRNA_DIR / "subclusters"

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
