# ============================================================================
# Purpose:      Figure 1F: SARS-CoV-2 spike IgG titers after the 3rd vaccine dose; SMM is split into Untreated and Treated columns. Same age- and sex-adjusted ANCOVA + JT framework as Figure 1B/C.
# Inputs:       data/elisa/elisa_spike_post3rd.csv.
# Outputs:      figures/Figure1F.png (and matching PDF).
# Dependencies: R + tidyverse, rstatix, ggpubr; sources ../config.R.
# ============================================================================
source("../config.R")
library(tidyverse)
library(rstatix)
library(ggpubr)

# Read ELISA titers after 3rd dose
Post_3rd_dose_df <- read.csv(
  file.path(DATA_DIR, "elisa", "elisa_spike_post3rd.csv")
)

df <- Post_3rd_dose_df %>%
  filter(Disease %in% c("Healthy", "MGUS", "IgM-MGUS", "SMM", "MM"))
df$Disease <- ifelse(df$Disease == "Healthy", "HD", as.character(df$Disease))
df$Disease_Recoded <- ifelse(df$Disease %in% c("IgM-MGUS"), "MGUS", as.character(df$Disease))

# Keep only first timepoint per individual
df <- df %>%
  group_by(Common_ID) %>%
  filter(row_number(ELISA_Titer) == 1) %>%
  ungroup()

# split SMM into SMM:Untreated and
# SMM:Treated using Ever_treated. MGUS uniformly treatment-naive. MM in IMPACT is by
# recruitment-cohort design previously or actively treated; not split.
df$Disease_Recoded <- ifelse(
  df$Disease_Recoded == "SMM" & df$Ever_treated == "Yes", "SMM:Treated",
  ifelse(df$Disease_Recoded == "SMM", "SMM:Untreated", as.character(df$Disease_Recoded))
)

# Wilcoxon tests vs HD with BH correction
df$Disease_Recoded <- factor(df$Disease_Recoded,
                             levels = c("HD", "MGUS", "SMM:Untreated", "SMM:Treated", "MM"))

# Plot
plotdf <- df
counts <- table(plotdf$Disease_Recoded)
# SMM:Untreated keeps tomato2; SMM:Treated darker red. Order matches factor levels.
color_palette <- c("steelblue", "orange", "tomato2", "#8B0000", "chartreuse4")

p <- ggplot(plotdf, aes(x = Disease_Recoded, y = ELISA_Titer, fill = Disease_Recoded)) +
  geom_violin(alpha = 0.5, position = position_dodge(width = .75), size = 1,
              show.legend = FALSE, color = NA) +
  geom_boxplot(lwd = 1, alpha = 0.5, outlier.size = -1, show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.2), shape = 21, alpha = 0.8,
             show.legend = FALSE) +
  scale_fill_manual(values = color_palette) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 12, family = FONT),
    axis.text.y = element_text(color = "black", size = 12, family = FONT),
    axis.title = element_text(size = 14, family = FONT),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    strip.text = element_text(size = 12, face = "plain", family = FONT),
    strip.background = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5, family = FONT)
  ) +
  xlab("") +
  ylab(bquote("ELISA Titer"["(OD450nm-570nm)"])) +
  ggtitle("After 3rd Dose") +
  scale_x_discrete(labels = function(x) {
    pretty <- gsub(":", "\n", x)
    paste0(pretty, "\n(n=", counts[x], ")")
  })

# Add q-value brackets with effect sizes
stat.test <- plotdf %>%
  rstatix::wilcox_test(ELISA_Titer ~ Disease_Recoded)
stat.test <- stat.test %>%
  filter(group1 == "HD") %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# Rank-biserial effect sizes for HD comparisons
eff <- plotdf %>%
  wilcox_effsize(ELISA_Titer ~ Disease_Recoded) %>%
  filter(group1 == "HD")

# AGE+SEX-adjusted q (rank-based ANCOVA), BH across HD contrasts.
# 4 HD-vs-each-group contrasts (MGUS, SMM:Untreated, SMM:Treated, MM).
.dd <- plotdf; .dd$rk <- rank(.dd$ELISA_Titer)
.cf <- coef(summary(lm(rk ~ Disease_Recoded + Sex + Age, data = .dd)))
.adj_q <- p.adjust(c(MGUS              = .cf["Disease_RecodedMGUS","Pr(>|t|)"],
                     `SMM:Untreated`   = .cf["Disease_RecodedSMM:Untreated","Pr(>|t|)"],
                     `SMM:Treated`     = .cf["Disease_RecodedSMM:Treated","Pr(>|t|)"],
                     MM                = .cf["Disease_RecodedMM","Pr(>|t|)"]), method = "BH")

stat.test <- stat.test %>%
  add_xy_position(x = "Disease_Recoded") %>%
  left_join(eff %>% select(group1, group2, effsize), by = c("group1", "group2")) %>%
  mutate(p.adj = .adj_q[as.character(group2)]) %>%
  mutate(
    q_fmt = ifelse(p.adj < 0.01,
                   paste0("q=", formatC(p.adj, format = "e", digits = 1)),
                   paste0("q=", round(p.adj, 2))),
    bracket_label = ifelse(p.adj < 0.1,
                           paste0(q_fmt, ", r=", round(abs(effsize), 2)),
                           q_fmt))
# tight bracket layout - see Figure1B.R for rationale
.ymax <- max(plotdf$ELISA_Titer, na.rm = TRUE)
.bracket_order <- c("MGUS" = 1, "SMM:Untreated" = 2, "SMM:Treated" = 3, "MM" = 4)
.step <- 0.42
stat.test$y.position <- .ymax + 0.25 +
  (.bracket_order[as.character(stat.test$group2)] - 1) * .step

# age- and sex-adjusted Jonckheere-Terpstra ordered-trend test.
.jt_df <- plotdf[plotdf$Disease_Recoded %in% c("HD", "MGUS", "SMM:Untreated"), ]
.jt <- jt_test_residuals_age_sex(.jt_df, "ELISA_Titer", "Disease_Recoded",
                                 order = c("HD", "MGUS", "SMM:Untreated"))
.jt_y    <- .ymax + 0.25 + 4 * .step + 0.40
.ylim_hi <- .jt_y + 0.55

bxp_padj <- p +
  stat_pvalue_manual(family = FONT, stat.test, label = "bracket_label",
                     tip.length = 0.02, size = 3.6,
                     inherit.aes = FALSE, hide.ns = FALSE) +
  annotate("segment", x = 1, xend = 3, y = .jt_y, yend = .jt_y, linewidth = 0.6) +
  annotate("segment", x = 1, xend = 1, y = .jt_y - 0.06, yend = .jt_y, linewidth = 0.6) +
  annotate("segment", x = 3, xend = 3, y = .jt_y - 0.06, yend = .jt_y, linewidth = 0.6) +
  annotate("text", family = FONT, x = 2.5, y = .jt_y + 0.30, label = jt_label(.jt),
           size = 3.8, fontface = "italic") +
  coord_cartesian(ylim = c(NA, .ylim_hi), clip = "off")

save_figure(bxp_padj, "Figure1F", width = 5.0, height = 4.0)
