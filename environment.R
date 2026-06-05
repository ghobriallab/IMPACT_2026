# ============================================================================
# Purpose:      One-shot R package installer for the figure scripts. Installs the union of CRAN packages used across Figure1-5 and SupplementalFigures.
# Inputs:       None (idempotent: skips packages already available).
# Outputs:      Installed R packages: tidyverse, ggplot2, ggpubr, rstatix, cowplot, ggrepel, plotrix, RColorBrewer, janitor, lme4, lmerTest, ineq, survival.
# Dependencies: R 4.0+ with internet access to a CRAN mirror.
# ============================================================================
# Install R packages required for IMPACT figure scripts

packages <- c(
  "tidyverse", "ggplot2", "ggpubr", "rstatix", "cowplot",
  "ggrepel", "plotrix", "RColorBrewer", "janitor",
  "lme4", "lmerTest", "ineq", "survival"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

cat("All packages installed.\n")
