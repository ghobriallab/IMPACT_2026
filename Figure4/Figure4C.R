source("../config.R")
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(plotrix)
library(rstatix)

# Load ClusTCR clonotype proportion data
clonotype_prop_df <- readRDS(file.path(DATA_DIR, "tcr", "tcr_clonotype_proportions.rds"))
average_prop <- clonotype_prop_df %>%
  group_by(PatientID, VaccineTimepoint, DiseaseStatus) %>%
  summarise_all(list(mean = "mean", sd = "sd"))

# Extract COVID/CEF columns
average_prop_covid <- average_prop[, c(1:7, 11:14)]

# Pre-vaccine
average_prop_pre <- average_prop_covid %>% filter(VaccineTimepoint == 1)
average_prop_pre <- average_prop_pre[, c(1, 2, 3, 4, 6, 8, 10)]
colnames(average_prop_pre)[4:7] <- c("CEF_mean_prop", "COVID_mean_prop", "CEF_SD_prop", "COVID_SD_prop")

# Post-vaccine
average_prop_post <- average_prop_covid %>% filter(VaccineTimepoint == 2)
average_prop_post <- average_prop_post[, c(1, 2, 3, 5, 7, 9, 11)]
colnames(average_prop_post)[4:7] <- c("CEF_mean_prop", "COVID_mean_prop", "CEF_SD_prop", "COVID_SD_prop")

# Combine and prepare for plotting
plot_df <- rbind(average_prop_pre, average_prop_post)
plot_df$VaccineTimepoint <- factor(plot_df$VaccineTimepoint, levels = c("1", "2"), ordered = TRUE)
plot_df$dummy_timepoint <- as.numeric(plot_df$VaccineTimepoint)

set.seed(123)
plot_df$add_noise <- plot_df$dummy_timepoint + runif(nrow(plot_df), -0.15, 0.15)

# Filter to paired samples only
plot_paired_df <- plot_df %>%
  group_by(PatientID) %>%
  filter(any(VaccineTimepoint == "1"), any(VaccineTimepoint == "2")) %>%
  ungroup()

# Convert to percentages
plot_paired_df$DiseaseStatus <- factor(plot_paired_df$DiseaseStatus, levels = c("HD", "MGUS", "SMM:Untreated"),
                                       labels = c("HD", "MGUS", "SMM"))
plot_paired_df <- plot_paired_df %>% mutate_at(vars(4:7), ~ . * 100)

# Disease labels with sample sizes
counts <- as.data.frame(table(plot_paired_df$VaccineTimepoint, plot_paired_df$DiseaseStatus))
counts <- counts[!duplicated(counts$Var2), ]
disease_names <- c("HD" = paste0("HD\n (n=", counts$Freq[counts$Var2 == "HD"], ")"),
                   "MGUS" = paste0("MGUS\n (n=", counts$Freq[counts$Var2 == "MGUS"], ")"),
                   "SMM" = paste0("SMM\n (n=", counts$Freq[counts$Var2 == "SMM"], ")"))

# CEF-specific TCR clonotypes
p_cef <- ggplot(plot_paired_df, aes(x = VaccineTimepoint, y = CEF_mean_prop, fill = VaccineTimepoint)) +
  geom_violin(alpha = 0.5, position = position_dodge(width = .75), linewidth = 1, show.legend = FALSE, color = NA) +
  geom_boxplot(alpha = 0.5, show.legend = FALSE) +
  geom_line(aes(x = add_noise, y = CEF_mean_prop, group = PatientID),
            color = "black", linewidth = 0.3, alpha = 0.5, show.legend = FALSE) +
  geom_pointrange(aes(x = add_noise, ymin = CEF_mean_prop - CEF_SD_prop, ymax = CEF_mean_prop + CEF_SD_prop),
                  shape = 21, size = 0.5, alpha = 1, show.legend = FALSE, linetype = "dotted") +
  scale_fill_manual(values = c("steelblue", "tomato2")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill = NA, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 13, color = "black"),
        strip.text = element_text(size = 15, face = "plain"),
        strip.background = element_rect(fill = NA)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  scale_x_discrete(breaks = c(1, 2), labels = c("Pre-Vx", "Post-Vx")) +
  facet_wrap(DiseaseStatus ~ ., labeller = as_labeller(disease_names)) +
  ylab("CEF-specific clonotypes (%)") +
  xlab("")

# Paired Wilcoxon test per disease group
stat_cef <- plot_paired_df %>%
  group_by(DiseaseStatus) %>%
  wilcox_test(CEF_mean_prop ~ VaccineTimepoint, paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance() %>%
  add_xy_position(x = "VaccineTimepoint")
stat_cef$p <- signif(stat_cef$p, 2)
stat_cef$p_label <- paste0("p = ", stat_cef$p)

p_cef_final <- p_cef +
  stat_pvalue_manual(stat_cef, label = "p_label", tip.length = 0.02, bracket.nudge.y = 0.01,
                     inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

ggsave(file.path(FIGURES_DIR, "Figure4C.png"), plot = p_cef_final, dpi = 300, units = "in", width = 6, height = 3.15)
