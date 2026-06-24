# ============================================================================
# Purpose:      Shared configuration for R figure scripts: data, repo and figure-output paths, color palette, and the age+sex-adjusted Jonckheere-Terpstra helper.
# Inputs:       Environment variable expectations; the user edits SCRNA_DIR to the directory containing scRNAseq_IMPACT_Zenodo.h5ad.
# Outputs:      DATA_DIR, FIGURES_DIR, COLORS, jt_test_residuals_age_sex() (sourced by every R figure script).
# Dependencies: R base; tryCatch; rstudioapi (optional, command-line Rscript falls back to getwd()).
# ============================================================================
# IMPACT — Configuration
# Edit these paths to match your local environment.

# Base directory containing scRNA-seq h5ad files (from GEO/Zenodo)
SCRNA_DIR <- "/path/to/scrnaseq_data"

# Base directory for this repository (auto-detected; works in RStudio AND command-line Rscript)
REPO_DIR <- tryCatch(
  dirname(rstudioapi::getSourceEditorContext()$path),
  error = function(e) ""
)
if (is.null(REPO_DIR) || REPO_DIR == "") {
  REPO_DIR <- getwd()
}

# Data files (download from Zenodo: https://doi.org/10.5281/zenodo.18989222)
DATA_DIR <- file.path(REPO_DIR, "..", "data")

# Output directory for figures
FIGURES_DIR <- file.path(REPO_DIR, "..", "figures")
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)

# age- and sex-adjusted Jonckheere-Terpstra
# ordered-trend test on residuals of rank(outcome) ~ Age + Sex.
# Returns list(z, p, n_per_group, n_used).  No external package needed.
jt_test_residuals_age_sex <- function(df, outcome, disease_col, order,
                                       age_col = "Age", sex_col = "Sex") {
  d <- df[complete.cases(df[, c(outcome, disease_col, age_col, sex_col)]), , drop = FALSE]
  d <- d[d[[disease_col]] %in% order, , drop = FALSE]
  if (nrow(d) < 3) return(list(z = NA_real_, p = NA_real_, nper = NULL, n = nrow(d)))
  d$.rk <- rank(d[[outcome]])
  # Residuals of rank ~ Age + Sex; categorical Sex via factor
  fmla <- as.formula(paste0(".rk ~ ", age_col, " + factor(", sex_col, ")"))
  d$.resid <- residuals(lm(fmla, data = d))
  groups <- lapply(order, function(g) d$.resid[d[[disease_col]] == g])
  groups <- groups[lengths(groups) > 0]
  k <- length(groups); N <- sum(lengths(groups))
  U <- 0
  for (i in seq_len(k - 1)) for (j in (i + 1):k) {
    xi <- groups[[i]]; xj <- groups[[j]]
    U <- U + sum(outer(xj, xi, `>`)) + 0.5 * sum(outer(xj, xi, `==`))
  }
  nis <- lengths(groups)
  EU <- (N^2 - sum(nis^2)) / 4
  VU <- (N^2 * (2*N + 3) - sum(nis^2 * (2*nis + 3))) / 72
  z <- (U - EU) / sqrt(VU)
  p <- 2 * (1 - pnorm(abs(z)))
  list(z = z, p = p,
       nper = setNames(as.integer(nis), order[seq_along(nis)]),
       n = N)
}

# Helper to format the JT label for on-panel annotation.
jt_label <- function(jt) {
  if (is.na(jt$z) || is.na(jt$p)) return(NA_character_)
  p_str <- if (jt$p < 0.001) sprintf("%.1e", jt$p) else sprintf("%.3f", jt$p)
  sprintf("Jonckheere-Terpstra: z=%.2f, p=%s", jt$z, p_str)
}

# Color palette (consistent across all figures)
COLORS <- c(
  "HD" = "#4DBBD5",
  "MGUS" = "#F39B7F",
  "SMM" = "#E64B35",
  "SMM (Untreated)" = "#E64B35",
  "SMM (Treated)" = "#91D1C2",
  "MM" = "#3C5488"
)
