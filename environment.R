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
