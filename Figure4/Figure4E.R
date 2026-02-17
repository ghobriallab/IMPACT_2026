source("../config.R")
library(dplyr)
library(tidyr)
library(rstatix)
library(ggpubr)

# Read ELISPOT data
raw_elispot_data <- read.csv(file.path(DATA_DIR, "elisa", "elispot_spike_cef.csv"))
raw_elispot_data <- raw_elispot_data[, c(1:11, 18:20)]

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

color_palette <- c("steelblue", "tomato2")

p <- ggplot(long_df, aes(x = Label, y = mean, fill = Disease)) +
  geom_violin(aes(Label, mean, fill = Disease),
              alpha = 0.5, position = position_dodge(width = .75),
              size = 1, show.legend = FALSE, color = NA) +
  geom_boxplot(alpha = 0.5, outlier.size = -1, show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.2), shape = 21, size = 2.5, show.legend = FALSE) +
  scale_fill_manual(values = color_palette) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 13),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title = element_text(size = 18),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        strip.text = element_text(size = 18, face = "plain"),
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
stat.test$p_label <- paste0("p = ", signif(stat.test$p, 2))

final_plot <- p +
  stat_pvalue_manual(stat.test, label = "{p_label}",
                     tip.length = 0.02, bracket.nudge.y = c(0.1),
                     size = 5.5, inherit.aes = FALSE, hide.ns = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)), n.breaks = 8)

ggsave(file.path(FIGURES_DIR, "Figure4E.png"), plot = final_plot, dpi = 300, units = "in", height = 4.5, width = 5)
