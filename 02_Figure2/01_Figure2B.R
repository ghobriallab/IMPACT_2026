# ============================================================================
# Purpose:      Figure 2B: MMR (Measles, Mumps, Rubella) antibody titers in HD vs SMM, adjusted for age (Age_Range midpoint) and sex via rank-based ANCOVA.
# Inputs:       data/elisa/elisa_mmr.csv (de-identified MMR ELISA).
# Outputs:      figures/Figure2B.png (and matching PDF and SVG).
# Dependencies: R + tidyverse, rstatix, ggpubr; sources ../config.R for FONT and save_figure().
# ============================================================================
source("../config.R")
library(tidyverse)
library(rstatix)
library(ggpubr)

# Read MMR ELISA data
MMR_df <- read.csv(file.path(DATA_DIR, "elisa", "elisa_mmr.csv"))

colnames(MMR_df)[2] <- "DiseaseStatus"
colnames(MMR_df)[4] <- "Rubella"
colnames(MMR_df)[6] <- "Mumps"
colnames(MMR_df)[8] <- "Measles"

# Reshape to long format
MMR_df_long <- gather(MMR_df, Infection, Titer_Measurement, 1:14,
                      -DiseaseStatus, -Age_Range, -Sample_ID, -Sex, -Race,
                      -FLC_ratio, -M_spike, -BM_PCs,
                      -Rubella_pos_neg, -Mumps_pos_neg, -Measles_pos_neg,
                      factor_key = TRUE)
MMR_df_long$Infection <- factor(MMR_df_long$Infection, levels = c("Measles", "Mumps", "Rubella"))

# --- age (Age_Range midpoint) + sex covariate adjustment ---
MMR_df_long$DiseaseStatus <- factor(MMR_df_long$DiseaseStatus, levels = c("HD", "SMM"))
.age_mid <- function(x) {
  m <- c("30-39"=34.5, "40-49"=44.5, "50-59"=54.5, "60-69"=64.5, "70-79"=74.5, "80-89"=84.5)
  unname(m[as.character(x)])
}
MMR_df_long$Age_mid <- .age_mid(MMR_df_long$Age_Range)

color_palette <- c("steelblue", "tomato2")

# Compute per-group sample counts (divide by 3 because long format has 3 infections per sample)
counts <- table(MMR_df_long$DiseaseStatus)
counts_numeric <- as.numeric(counts / 3)
names(counts_numeric) <- names(counts)
counts_numeric <- as.table(counts_numeric)

# Wilcoxon tests per infection, BH-corrected
stat.test <- MMR_df_long %>%
  group_by(Infection) %>%
  rstatix::wilcox_test(Titer_Measurement ~ DiseaseStatus) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>%
  add_xy_position(x = "DiseaseStatus", scales = "free")

# Compute y-positions from initial plot build
p_tmp <- ggplot(MMR_df_long, aes(x = DiseaseStatus, y = Titer_Measurement, fill = DiseaseStatus)) +
  geom_violin(alpha = 0.5, size = 1, show.legend = FALSE, color = NA) +
  scale_y_log10() +
  facet_wrap(~Infection, scales = "free")
ypos_pval <- c(ggplot_build(p_tmp)$layout$panel_scales_y[[1]]$range$range[2],
               ggplot_build(p_tmp)$layout$panel_scales_y[[2]]$range$range[2],
               ggplot_build(p_tmp)$layout$panel_scales_y[[3]]$range$range[2])
stat.test$y.position <- ypos_pval

stat.test$p.adj.signif <- ifelse(stat.test$p.adj.signif == "ns", "", stat.test$p.adj.signif)

# replace displayed q with AGE+SEX-adjusted q (rank-based ANCOVA), BH across infections
.adj_p <- sapply(levels(MMR_df_long$Infection), function(inf) {
  s <- subset(MMR_df_long, Infection == inf); s$rk <- rank(s$Titer_Measurement)
  coef(summary(lm(rk ~ DiseaseStatus + Age_mid + Sex, data = s)))["DiseaseStatusSMM", "Pr(>|t|)"]
})
.adj_q <- p.adjust(.adj_p, method = "BH")
stat.test$p.adj <- .adj_q[as.character(stat.test$Infection)]
stat.test$p.adj <- paste0("q=", sprintf("%.2g", stat.test$p.adj))

# Final plot
p <- ggplot(MMR_df_long, aes(x = DiseaseStatus, y = Titer_Measurement, fill = DiseaseStatus)) +
  geom_violin(alpha = 0.5, position = position_dodge(width = .75), size = 1,
              show.legend = FALSE, color = NA) +
  geom_boxplot(lwd = 1, alpha = 0.5, outlier.size = -1, show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.7, shape = 21,
             show.legend = FALSE) +
  scale_fill_manual(values = color_palette) +
  theme(
    text = element_text(family = FONT),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 12, family = FONT),
    axis.text.y = element_text(color = "black", size = 12, family = FONT),
    axis.title = element_text(size = 14, family = FONT),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    strip.text = element_text(size = 14, face = "plain", family = FONT),
    strip.background = element_blank()
  ) +
  xlab("") +
  ylab("Antibody Titer (IU/ml)") +
  scale_y_log10(expand = expansion(mult = c(0.05, 0.15))) +
  stat_pvalue_manual(stat.test, label = "{p.adj}", family = FONT,
                     tip.length = 0.02, bracket.nudge.y = c(0.2),
                     inherit.aes = FALSE, hide.ns = FALSE) +
  facet_wrap(~Infection, scales = "free") +
  scale_x_discrete(labels = function(x) paste0(x, "\n(n=", counts_numeric[x], ")"))

save_figure(p, "Figure2B", width = 7, height = 3)
