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
df <- df[order(as.Date(df$Date, format = "%Y/%m/%d")), ]

# Keep only first timepoint per individual
df <- df %>%
  group_by(Common_ID) %>%
  filter(row_number(ELISA_Titer) == 1) %>%
  ungroup()

# Wilcoxon tests vs HD with BH correction
df$Disease_Recoded <- factor(df$Disease_Recoded, levels = c("HD", "MGUS", "SMM", "MM"))
stat.test_bundled <- df %>%
  rstatix::wilcox_test(ELISA_Titer ~ Disease_Recoded)
stat.test_bundled <- stat.test_bundled %>%
  filter(group1 == "HD") %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# Plot
color_palette <- c("steelblue", "orange", "tomato2", "chartreuse4", "deeppink")
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
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title = element_text(size = 12),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    strip.text = element_text(size = 12, face = "plain"),
    strip.background = element_blank(),
    plot.title = element_text(size = 12)
  ) +
  xlab("") +
  ylab(bquote(ELISA_Titer["(OD450nm-570nm)"])) +
  ggtitle("2 weeks - 2 months from 2nd dose") +
  scale_x_discrete(labels = function(x) paste0(x, "\n(n=", counts[x], ")"))

# Rank-biserial effect sizes for HD comparisons
eff <- df %>%
  wilcox_effsize(ELISA_Titer ~ Disease_Recoded) %>%
  filter(group1 == "HD")

# Add q-value brackets (HD vs each group), with effect size when q < 0.1
stat.test_bundled <- stat.test_bundled %>%
  add_xy_position(x = "Disease_Recoded", step.increase = c(0.5))
stat.test_bundled_comparison_with_HD <- stat.test_bundled %>%
  filter(group1 == "HD") %>%
  left_join(eff %>% select(group1, group2, effsize), by = c("group1", "group2")) %>%
  mutate(
    q_fmt = ifelse(p.adj < 0.01,
                   paste0("q=", formatC(p.adj, format = "e", digits = 1)),
                   paste0("q=", round(p.adj, 2))),
    bracket_label = ifelse(p.adj < 0.1,
                           paste0(q_fmt, ", r=", round(abs(effsize), 2)),
                           q_fmt))

bxp_padj <- p +
  stat_pvalue_manual(stat.test_bundled_comparison_with_HD, label = "bracket_label",
                     tip.length = 0.02, bracket.nudge.y = c(0.2),
                     inherit.aes = FALSE, hide.ns = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

ggsave(file.path(FIGURES_DIR, "Figure1B.png"), bxp_padj, dpi = 300, unit = "in", width = 4, height = 3.5)
