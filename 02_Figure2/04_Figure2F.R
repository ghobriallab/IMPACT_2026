# ============================================================================
# Purpose:      Figure 2F: post-vaccination APRIL (TNFSF13) vs serum M-spike Spearman correlation in treatment-naive SMM patients with a quantifiable M-spike. Includes a nonparametric bootstrap CI on rho and a leave-one-out sensitivity analysis.
# Inputs:       data/olink/olink_summary_tumor_burden.csv.
# Outputs:      figures/Figure2F.png (and matching PDF and SVG) with the Spearman rho, p-value and bootstrap CI annotated.
# Dependencies: R + tidyverse; sources ../config.R for FONT and save_figure().
# ============================================================================
source("../config.R")
library(tidyverse)

# Read Olink data with tumor burden (M-spike)
bm <- read.csv(file.path(DATA_DIR, "olink", "olink_summary_tumor_burden.csv"),
               check.names = FALSE)

# Filter: Post-vaccination, TNFSF13 (APRIL), SMM patients, QC pass
bm_post <- bm[bm$`Timepoint (red_within 14 days, orange_with pre-vax sample)` == "After_Vax", ]
bm_post <- bm_post[bm_post$Cytokine == "TNFSF13" &
                    bm_post$Diagnosis_when_fully_vaccinated %in% c("SMM") &
                    bm_post$QC_1_pass_0_warning == 1, ]
# require a recorded quantitative M-spike for inclusion in the
# correlation; this is why the plotted n=10 rather than the 18 SMM patients in the Olink
# cohort (M-spike was not quantifiable in the remaining 8).
bm_post <- bm_post[!is.na(bm_post$M_spike) & !is.na(bm_post$Level), ]

# Spearman correlation (robust to non-normal distributions)
spearman_test <- cor.test(bm_post$M_spike, bm_post$Level, method = "spearman", exact = FALSE)
rho_value <- signif(spearman_test$estimate, 2)
p_value <- signif(spearman_test$p.value, 1)
n_used <- nrow(bm_post)

# robustness diagnostics. Nonparametric bootstrap
# (10,000 resamples) for a 95% CI on the Spearman rho, and a leave-one-out sensitivity analysis.
# Numbers also reported in the legend and Methods. Seed fixed for reproducibility.
set.seed(2026)
B <- 10000
boot_rho <- replicate(B, {
  idx <- sample.int(n_used, n_used, replace = TRUE)
  if (length(unique(idx)) < 3) return(NA_real_)
  cor(bm_post$M_spike[idx], bm_post$Level[idx], method = "spearman")
})
boot_rho <- boot_rho[!is.na(boot_rho)]
# Format CI as two-decimal strings so "-1.00 to -0.70" prints with trailing zeros (signif() drops them).
ci_lo <- sprintf("%.2f", quantile(boot_rho, 0.025))
ci_hi <- sprintf("%.2f", quantile(boot_rho, 0.975))
rho_str <- sprintf("%.2f", spearman_test$estimate)

p <- ggplot(bm_post) +
  geom_point(aes(M_spike, Level), alpha = 0.75, size = 3, shape = 21,
             color = "black", fill = "steelblue") +
  geom_smooth(aes(M_spike, Level), method = "lm", color = "tomato2",
              fill = "lightblue", level = 0.95) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    text = element_text(family = FONT),
    axis.text = element_text(size = 12, color = "black", family = FONT),
    axis.title = element_text(size = 14, family = FONT),
    plot.title = element_text(size = 14, hjust = 0.5, family = FONT)
  ) +
  xlab("M-spike (g/dL)") +
  # Parentheses inside a plotmath GROUP get typeset from a symbol font, which breaks the all-Arial
  # rule. Quoting them keeps them in Arial while TNFSF13 stays italic per the gene-name convention.
  ylab(expression("Post-Vx " * italic("TNFSF13") * " (APRIL) level")) +
  ggtitle(paste0("SMM patients (n=", n_used, ")")) +
  # On-figure annotation: Spearman rho, bootstrap 95% CI on its own line, then p on the third line.
  # n is shown only in the title (no need to repeat). Right-justified upper-right.
  annotate("text", x = Inf, y = Inf,
           label = paste0("Spearman ρ = ", rho_str, "\n",
                          "95% CI ", ci_lo, " to ", ci_hi, "\n",
                          "p = ", p_value),
           color = "black", size = 3.4, hjust = 1.05, vjust = 1.3, family = FONT,
           lineheight = 1.05) +
  # Expand y-axis upper limit so the annotation sits above all data points; do NOT occlude the
  # reviewer-flagged influential point at M-spike=0 (highest APRIL).
  scale_y_continuous(limits = c(NA, 42000))

save_figure(p, "Figure2F", width = 3.2, height = 4.8)
