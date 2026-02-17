# IMPACT 2026 — Configuration
# Edit these paths to match your local environment.

# Base directory containing scRNA-seq h5ad files (from GEO/Zenodo)
SCRNA_DIR <- "/path/to/scrnaseq_data"

# Base directory for this repository (auto-detected)
REPO_DIR <- dirname(rstudioapi::getSourceEditorContext()$path)
if (is.null(REPO_DIR) || REPO_DIR == "") {
  REPO_DIR <- getwd()
}

# Data files (download from Zenodo: https://doi.org/10.5281/zenodo.18989223)
DATA_DIR <- file.path(REPO_DIR, "..", "data")

# Output directory for figures
FIGURES_DIR <- file.path(REPO_DIR, "..", "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

# Color palette (consistent across all figures)
COLORS <- c(
  "HD" = "#4DBBD5",
  "MGUS" = "#F39B7F",
  "SMM" = "#E64B35",
  "SMM (Untreated)" = "#E64B35",
  "SMM (Treated)" = "#91D1C2",
  "MM" = "#3C5488"
)
