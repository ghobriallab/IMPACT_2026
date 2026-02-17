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

# Wilcoxon tests vs HD with BH correction
df$Disease_Recoded <- factor(df$Disease_Recoded, levels = c("HD", "MGUS", "SMM", "MM"))

# Plot
plotdf <- df
counts <- table(plotdf$Disease_Recoded)
color_palette <- c("steelblue", "orange", "tomato2", "chartreuse4", "deeppink")

p <- ggplot(plotdf, aes(x = Disease_Recoded, y = ELISA_Titer, fill = Disease_Recoded)) +
  geom_violin(alpha = 0.5, position = position_dodge(width = .75), size = 1,
              show.legend = FALSE, color = NA) +
  geom_boxplot(lwd = 1, alpha = 0.5, outlier.size = -1, show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.2), shape = 21, alpha = 0.8,
             show.legend = FALSE) +
  scale_fill_manual(values = color_palette) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title = element_text(size = 14),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    strip.text = element_text(size = 12, face = "plain"),
    strip.background = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5)
  ) +
  xlab("") +
  ylab(bquote(ELISA ~ Titer["(OD450nm-570nm)"])) +
  ggtitle("After 3rd Dose") +
  scale_x_discrete(labels = function(x) paste0(x, "\n(n=", counts[x], ")"))

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

stat.test <- stat.test %>%
  add_xy_position(x = "Disease_Recoded") %>%
  left_join(eff %>% select(group1, group2, effsize), by = c("group1", "group2")) %>%
  mutate(
    q_fmt = ifelse(p.adj < 0.01,
                   paste0("q=", formatC(p.adj, format = "e", digits = 1)),
                   paste0("q=", round(p.adj, 2))),
    bracket_label = ifelse(p.adj < 0.1,
                           paste0(q_fmt, ", r=", round(abs(effsize), 2)),
                           q_fmt))
stat.test$y.position <- c(3.2, 3.5, 3.8)

bxp_padj <- p +
  stat_pvalue_manual(stat.test, label = "bracket_label",
                     tip.length = 0.02, size = 3,
                     inherit.aes = FALSE, hide.ns = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

ggsave(file.path(FIGURES_DIR, "Figure1F.png"), bxp_padj, dpi = 300, unit = "in", width = 4, height = 4)
