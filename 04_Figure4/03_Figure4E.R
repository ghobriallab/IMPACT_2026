# ============================================================================
# Purpose:      Figure 4E: IFN-gamma ELISPOT for SARS-CoV-2 spike (SARS) and CERI recall pools, in HD vs SMM, technical-duplicate-averaged and DMSO-normalized; age- and sex-adjusted rank-based ANCOVA.
# Inputs:       data/elisa/elispot_spike_cef.csv (de-identified ELISPOT well-level counts).
# Outputs:      figures/Figure4E.png (and matching PDF and SVG).
# Dependencies: R + dplyr, tidyr, rstatix, ggpubr; sources ../config.R for FONT and save_figure().
# ============================================================================
source("../config.R")
library(dplyr)
library(tidyr)
library(rstatix)
library(ggpubr)

# Read ELISPOT data
raw_elispot_data <- read.csv(file.path(DATA_DIR, "elisa", "elispot_spike_cef.csv"))
# select needed columns by name (positional index c(1:11,18:20) referenced a
# non-existent column 20 in the shared 19-column file and errored).
raw_elispot_data <- raw_elispot_data[, c("Sample_ID", "Disease", "Age_Range", "Sex",
                                         "DMSO_1", "CERI_1", "SARS_1", "DMSO_2", "CERI_2", "SARS_2",
                                         "FLC_ratio", "M_spike", "BM_PCs")]

# Calculate ratios to DMSO
raw_elispot_data$Ratio_CERI_1 <- raw_elispot_data$CERI_1 / raw_elispot_data$DMSO_1
raw_elispot_data$Ratio_CERI_2 <- raw_elispot_data$CERI_2 / raw_elispot_data$DMSO_2
raw_elispot_data$Ratio_SARS_1 <- raw_elispot_data$SARS_1 / raw_elispot_data$DMSO_1
raw_elispot_data$Ratio_SARS_2 <- raw_elispot_data$SARS_2 / raw_elispot_data$DMSO_2

raw_elispot_data$mean_CERI_ratio <- rowMeans(raw_elispot_data[, c("Ratio_CERI_1", "Ratio_CERI_2")], na.rm = TRUE)
raw_elispot_data$mean_SARS_ratio <- rowMeans(raw_elispot_data[, c("Ratio_SARS_1", "Ratio_SARS_2")], na.rm = TRUE)

raw_elispot_data <- raw_elispot_data %>%
  group_by(Sample_ID) %>%
  mutate(SD_CERI_ratio = sd(unlist(select(cur_data(), Ratio_CERI_1:Ratio_CERI_2))))
raw_elispot_data <- raw_elispot_data %>%
  group_by(Sample_ID) %>%
  mutate(SD_SARS_ratio = sd(unlist(select(cur_data(), Ratio_SARS_1:Ratio_SARS_2))))

# Reshape to long format
df <- raw_elispot_data[, c("Sample_ID", "Disease", "mean_CERI_ratio", "mean_SARS_ratio", "SD_CERI_ratio", "SD_SARS_ratio")]
long_df <- df %>%
  pivot_longer(
    cols = starts_with(c("mean", "SD")),
    names_to = c(".value", "category"),
    names_pattern = "(mean|SD)_(.*)"
  )

# Add sample count labels
count_data <- long_df %>%
  group_by(category, Disease) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Label = paste0(Disease, "\n(n=", Count, ")"))

long_df <- long_df %>%
  left_join(count_data, by = c("category", "Disease"))

long_df$Label <- factor(long_df$Label, levels = c("HD\n(n=10)", "SMM\n(n=10)"))
long_df$category <- factor(long_df$category, levels = c("SARS_ratio", "CERI_ratio"))

# --- merge age (Age_Range midpoint) + sex for covariate adjustment ---
.age_mid <- function(x) {
  m <- c("30-39"=34.5, "40-49"=44.5, "50-59"=54.5, "60-69"=64.5, "70-79"=74.5, "80-89"=84.5)
  unname(m[as.character(x)])
}
.meta_as <- raw_elispot_data %>% group_by(Sample_ID) %>%
  summarise(Age_Range = first(Age_Range), Sex = first(Sex), .groups = "drop")
long_df <- left_join(long_df, .meta_as, by = "Sample_ID")
long_df$Age_mid <- .age_mid(long_df$Age_Range)
long_df$Disease <- factor(long_df$Disease, levels = c("HD", "SMM"))

color_palette <- c("steelblue", "tomato2")

p <- ggplot(long_df, aes(x = Label, y = mean, fill = Disease)) +
  geom_violin(aes(Label, mean, fill = Disease),
              alpha = 0.5, position = position_dodge(width = .75),
              size = 1, show.legend = FALSE, color = NA) +
  geom_boxplot(alpha = 0.5, outlier.size = -1, show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.2), shape = 21, size = 2.5, show.legend = FALSE) +
  scale_fill_manual(values = color_palette) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 13, family = FONT),
        axis.text.y = element_text(color = "black", size = 16, family = FONT),
        axis.title = element_text(size = 18, family = FONT),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(size = 18, face = "plain", family = FONT),
        strip.background = element_blank()) +
  xlab("") +
  ylab(NULL) +
  facet_wrap(~category, scales = "free_y", strip.position = "left",
             labeller = as_labeller(c(SARS_ratio = "SARS CoV2/DMSO", CERI_ratio = "CERI/DMSO"))) +
  theme(strip.background = element_blank(),
        strip.placement = "outside")

# Wilcoxon test per category (two independent tests, raw p-values)
stat.test <- long_df %>%
  group_by(category) %>%
  rstatix::wilcox_test(mean ~ Disease) %>%
  add_significance()
stat.test <- stat.test %>% add_xy_position(x = "Disease", scales = "free")

# replace displayed p with AGE+SEX-adjusted p (rank-based ANCOVA), per category
.adj_p <- sapply(levels(long_df$category), function(cat) {
  s <- subset(long_df, category == cat & !is.na(Sex) & !is.na(Age_mid)); s$rk <- rank(s$mean)
  cf <- coef(summary(lm(rk ~ Disease + Age_mid + Sex, data = s)))
  cf[grep("^Disease", rownames(cf))[1], "Pr(>|t|)"]
})
stat.test$p <- .adj_p[as.character(stat.test$category)]
stat.test$p_label <- paste0("p = ", signif(stat.test$p, 2))

final_plot <- p +
  stat_pvalue_manual(family = FONT, stat.test, label = "{p_label}",
                     tip.length = 0.02, bracket.nudge.y = c(0.1),
                     size = 5.5, inherit.aes = FALSE, hide.ns = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)), n.breaks = 8)

save_figure(final_plot, "Figure4E", width = 5, height = 4.5)
