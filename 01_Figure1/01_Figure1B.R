# ============================================================================
# Purpose:      Figure 1B: SARS-CoV-2 spike IgG titers at the early post-2nd-dose timepoint (2-8 weeks). Compares HD, MGUS, SMM and MM with age- and sex-adjusted rank-based ANCOVA and Jonckheere-Terpstra ordered-trend test.
# Inputs:       data/elisa/elisa_spike_post2nd.csv (de-identified ELISA titers).
# Outputs:      figures/Figure1B.png (and matching PDF).
# Dependencies: R + tidyverse, rstatix, ggpubr; sources ../config.R.
# ============================================================================
source("../config.R")
library(tidyverse)
library(rstatix)
library(ggpubr)

# Read ELISA titers after 2nd dose
Post_2nd_dose_df <- read.csv(
  file.path(DATA_DIR, "elisa", "elisa_spike_post2nd.csv")
)

df <- Post_2nd_dose_df %>%
  filter(Disease %in% c("Healthy", "MGUS", "IgM-MGUS", "SMM", "MM"))
df$Disease <- ifelse(df$Disease == "Healthy", "HD", as.character(df$Disease))
df$Disease_Recoded <- ifelse(df$Disease %in% c("IgM-MGUS"), "MGUS", as.character(df$Disease))

# Filter: 2 weeks to 2 months post-2nd dose
df <- df %>% filter(Days_post2nd <= 60 & Days_post2nd > 13)
df <- df[order(df$Days_post2nd), ]  # order by days post-2nd (no Date column in shared data)

# Keep only first timepoint per individual
df <- df %>%
  group_by(Common_ID) %>%
  filter(row_number(ELISA_Titer) == 1) %>%
  ungroup()

# split SMM into SMM:Untreated and
# SMM:Treated using Ever_treated. MGUS is uniformly treatment-naive across the cohort (no Yes).
# MM in IMPACT is clinically/by-recruitment-design previously or actively treated; we do not
# attempt to split MM since Ever_treated within MM is essentially collinear with disease.
df$Disease_Recoded <- ifelse(
  df$Disease_Recoded == "SMM" & df$Ever_treated == "Yes", "SMM:Treated",
  ifelse(df$Disease_Recoded == "SMM", "SMM:Untreated", as.character(df$Disease_Recoded))
)

# Wilcoxon tests vs HD with BH correction
df$Disease_Recoded <- factor(df$Disease_Recoded,
                             levels = c("HD", "MGUS", "SMM:Untreated", "SMM:Treated", "MM"))
stat.test_bundled <- df %>%
  rstatix::wilcox_test(ELISA_Titer ~ Disease_Recoded)
stat.test_bundled <- stat.test_bundled %>%
  filter(group1 == "HD") %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# Plot
# SMM:Untreated keeps the canonical tomato2; SMM:Treated darker red to visually
# group as the SMM family while remaining distinct. Color order matches factor levels above.
color_palette <- c("steelblue", "orange", "tomato2", "#8B0000", "chartreuse4")
plotdf <- df
counts <- table(plotdf$Disease_Recoded)

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
  ggtitle("2 weeks - 2 months from 2nd dose") +
  scale_x_discrete(labels = function(x) {
    # Wrap SMM:Untreated / SMM:Treated across two lines so axis labels do not overlap
    pretty <- gsub(":", "\n", x)
    paste0(pretty, "\n(n=", counts[x], ")")
  })

# Rank-biserial effect sizes for HD comparisons
eff <- df %>%
  wilcox_effsize(ELISA_Titer ~ Disease_Recoded) %>%
  filter(group1 == "HD")

# AGE+SEX-adjusted q (rank-based ANCOVA), BH across HD contrasts.
# Age at sampling (capped at 90 per HIPAA Safe Harbor) joined from the deid data; complete-case for the model.
# now 4 HD-vs-each-group contrasts (MGUS, SMM:Untreated, SMM:Treated, MM).
.dd <- df; .dd$rk <- rank(.dd$ELISA_Titer)
.cf <- coef(summary(lm(rk ~ Disease_Recoded + Sex + Age, data = .dd)))
.adj_q <- p.adjust(c(MGUS              = .cf["Disease_RecodedMGUS","Pr(>|t|)"],
                     `SMM:Untreated`   = .cf["Disease_RecodedSMM:Untreated","Pr(>|t|)"],
                     `SMM:Treated`     = .cf["Disease_RecodedSMM:Treated","Pr(>|t|)"],
                     MM                = .cf["Disease_RecodedMM","Pr(>|t|)"]), method = "BH")

# Add q-value brackets (HD vs each group), with effect size when q < 0.1
stat.test_bundled <- stat.test_bundled %>%
  add_xy_position(x = "Disease_Recoded")
stat.test_bundled_comparison_with_HD <- stat.test_bundled %>%
  filter(group1 == "HD") %>%
  left_join(eff %>% select(group1, group2, effsize), by = c("group1", "group2")) %>%
  mutate(p.adj = .adj_q[as.character(group2)]) %>%
  mutate(
    q_fmt = ifelse(p.adj < 0.01,
                   paste0("q=", formatC(p.adj, format = "e", digits = 1)),
                   paste0("q=", round(p.adj, 2))),
    bracket_label = ifelse(p.adj < 0.1,
                           paste0(q_fmt, ", r=", round(abs(effsize), 2)),
                           q_fmt))
# tight bracket layout - pairwise brackets stacked closely, JT bracket on top.
# coord_cartesian provides explicit y-limit so the panel isn't padded with empty white space
# above the annotations (scale_y_continuous expansion would auto-extend by mult * range).
.ymax <- max(plotdf$ELISA_Titer, na.rm = TRUE)
.bracket_order <- c("MGUS" = 1, "SMM:Untreated" = 2, "SMM:Treated" = 3, "MM" = 4)
.step <- 0.42                                        # readable but not airy
stat.test_bundled_comparison_with_HD$y.position <-
  .ymax + 0.25 + (.bracket_order[as.character(stat.test_bundled_comparison_with_HD$group2)] - 1) * .step

# age- and sex-adjusted Jonckheere-Terpstra ordered-trend test for
# monotonic disease progression HD < MGUS < SMM:Untreated. Bracket sits above the 4 pairwise
# brackets and spans columns 1-3.
.jt_df <- df[df$Disease_Recoded %in% c("HD", "MGUS", "SMM:Untreated"), ]
.jt <- jt_test_residuals_age_sex(.jt_df, "ELISA_Titer", "Disease_Recoded",
                                 order = c("HD", "MGUS", "SMM:Untreated"))
.jt_y    <- .ymax + 0.25 + 4 * .step + 0.40    # moderate gap above highest pairwise (MM)
.ylim_hi <- .jt_y + 0.55                       # tight upper limit (no extra empty space)

bxp_padj <- p +
  stat_pvalue_manual(family = FONT, stat.test_bundled_comparison_with_HD, label = "bracket_label",
                     tip.length = 0.02, size = 3.6,
                     inherit.aes = FALSE, hide.ns = FALSE) +
  annotate("segment", x = 1, xend = 3, y = .jt_y, yend = .jt_y, linewidth = 0.6) +
  annotate("segment", x = 1, xend = 1, y = .jt_y - 0.06, yend = .jt_y, linewidth = 0.6) +
  annotate("segment", x = 3, xend = 3, y = .jt_y - 0.06, yend = .jt_y, linewidth = 0.6) +
  annotate("text", family = FONT, x = 2.5, y = .jt_y + 0.30, label = jt_label(.jt),
           size = 3.8, fontface = "italic") +
  coord_cartesian(ylim = c(NA, .ylim_hi), clip = "off")

save_figure(bxp_padj, "Figure1B", width = 5.0, height = 4.0)
