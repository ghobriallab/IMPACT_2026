source("../config.R")
library(tidyverse)
library(rstatix)
library(ggpubr)

# Read Olink protein data
olink <- read.csv(file.path(DATA_DIR, "olink", "olink_paired_prepost.csv"), check.names = FALSE)
olink$Diagnosis <- olink$Diagnosis_when_fully_vaccinated
olink_meta <- olink[!duplicated(olink$Patient_ID), ]

# Remove NA and IgM-MGUS
olink <- olink[!is.na(olink$Level), ]
olink <- olink[olink$Diagnosis != "IgM-MGUS", ]

# Collapse replicate estimates by mean
olink$P_C_T <- paste0(olink$Patient_ID, ":", olink$Cytokine, ";", olink$Timepoint)
olink_2 <- aggregate(olink$Level, list(olink$P_C_T), mean, na.rm = TRUE)
colnames(olink_2) <- c("P_C_T", "Level")
olink_2$Cytokine <- gsub("^.*:", "", olink_2$P_C_T)
olink_2$Cytokine <- gsub(";.*$", "", olink_2$Cytokine)
olink_2$Timepoint <- gsub("^.*;", "", olink_2$P_C_T)
olink_2$Patient_ID <- gsub(":.*$", "", olink_2$P_C_T)
olink_3 <- aggregate(olink$QC_1_pass_0_warning, list(olink$P_C_T), max, na.rm = TRUE)
colnames(olink_3) <- c("P_C_T", "QC_1_pass_0_warning")
olink_2$QC_1_pass_0_warning <- olink_3$QC_1_pass_0_warning[match(olink_2$P_C_T, olink_3$P_C_T)]
olink_2$P_C_T <- NULL
olink_2$Diagnosis <- olink_meta$Diagnosis[match(olink_2$Patient_ID, olink_meta$Patient_ID)]

# QC filter: keep only samples that pass QC and have no breakthrough infection
data <- olink_2 %>% filter(QC_1_pass_0_warning == "1")
data <- data[data$Diagnosis != "IgM-MGUS", ]

# Identify cytokines with sufficient data across all groups
check_data <- as.data.frame(table(data$Cytokine, data$Diagnosis, data$Timepoint))
colnames(check_data) <- c("Cytokine", "Diagnosis", "Timepoint", "Freq")
check_data_wide <- reshape(check_data, idvar = c("Cytokine", "Timepoint"),
                           timevar = "Diagnosis", direction = "wide")
cytokines_with_enough_data <- check_data_wide %>%
  group_by(Timepoint) %>%
  filter(Freq.Healthy > 1 & Freq.MGUS > 1 & Freq.SMM > 1) %>%
  pull(Cytokine)
cytokines_with_enough_data <- unique(cytokines_with_enough_data)

data$Timepoint <- ifelse(data$Timepoint == "After_Vax", "Post-Vx", "Pre-Vx")
data$Timepoint <- factor(data$Timepoint, levels = c("Pre-Vx", "Post-Vx"))

# BH-corrected Wilcoxon tests across all cytokines (HD vs others)
# Pre-Vx
stat.test_prevax_w_MGUS <- data %>%
  filter(Timepoint == "Pre-Vx") %>%
  filter(Cytokine %in% cytokines_with_enough_data) %>%
  group_by(Timepoint, Cytokine) %>%
  rstatix::wilcox_test(Level ~ Diagnosis) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test_prevax_w_MGUS <- stat.test_prevax_w_MGUS %>% filter(!group1 %in% c("MGUS"))
stat.test_prevax_w_MGUS <- stat.test_prevax_w_MGUS %>%
  add_xy_position(x = "Diagnosis", step.increase = c(0.1))
stat.test_prevax_w_MGUS$p.adj <- paste0("q=", signif(stat.test_prevax_w_MGUS$p.adj, 2))

# Post-Vx
stat.test_postvax_w_MGUS <- data %>%
  filter(Timepoint == "Post-Vx") %>%
  filter(Cytokine %in% cytokines_with_enough_data) %>%
  group_by(Timepoint, Cytokine) %>%
  rstatix::wilcox_test(Level ~ Diagnosis) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test_postvax_w_MGUS <- stat.test_postvax_w_MGUS %>% filter(!group1 %in% c("MGUS"))
stat.test_postvax_w_MGUS <- stat.test_postvax_w_MGUS %>%
  add_xy_position(x = "Diagnosis", step.increase = c(0.1))
stat.test_postvax_w_MGUS$p.adj <- paste0("q=", signif(stat.test_postvax_w_MGUS$p.adj, 2))

new_stat.test <- rbind(stat.test_prevax_w_MGUS, stat.test_postvax_w_MGUS)

# --- Plot TNFSF13 (APRIL) ---
all_cyto <- unique(as.character(new_stat.test$Cytokine))
i <- which(all_cyto == "TNFSF13")
tmp <- data %>% filter(Cytokine == all_cyto[i])

color_palette <- c("steelblue", "orange", "tomato2", "chartreuse4", "deeppink")

# Per-group counts for axis labels
count_data <- tmp %>%
  group_by(Timepoint, Diagnosis) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Label = paste0(Diagnosis, "\n(n=", Count, ")"))
tmp <- tmp %>% left_join(count_data, by = c("Timepoint", "Diagnosis"))

tmp$Diagnosis <- factor(tmp$Diagnosis, levels = c("Healthy", "MGUS", "SMM"))

p <- ggplot(tmp, aes(x = Label, y = Level, fill = Diagnosis)) +
  geom_violin(alpha = 0.5, position = position_dodge(width = .75), size = 1,
              show.legend = FALSE, color = NA, outlier.shape = NA) +
  geom_boxplot(lwd = 1, alpha = 0.5, outlier.shape = NA, show.legend = FALSE) +
  geom_point(position = position_jitter(width = 0.2), shape = 21,
             outlier.shape = NA, outlier.size = 1, show.legend = FALSE) +
  scale_fill_manual(values = color_palette) +
  facet_wrap(~ Timepoint, scales = "free") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title = element_text(size = 14),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    strip.text = element_text(size = 14, face = "plain"),
    strip.background = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5)
  ) +
  xlab("") +
  ylab("TNFSF13 (APRIL) Level (NPX)")

# q-values for TNFSF13 from the BH-corrected full analysis
stat.test_cytokine_1 <- new_stat.test %>%
  filter(Cytokine == "TNFSF13") %>%
  filter(group1 == "Healthy")

# Bracket positions
stat.test_cytokine_1$xmin <- ifelse(stat.test_cytokine_1$group2 == "MGUS", 1, 1)
stat.test_cytokine_1$xmax <- ifelse(stat.test_cytokine_1$group2 == "MGUS", 2, 3)
max_level <- max(tmp$Level, na.rm = TRUE)
stat.test_cytokine_1$y.position <- ifelse(stat.test_cytokine_1$group2 == "MGUS",
                                          max_level * 1.05, max_level * 1.20)

# Format q-value labels
stat.test_cytokine_1$q_label <- ifelse(
  as.numeric(gsub("q=", "", stat.test_cytokine_1$p.adj)) < 0.01,
  paste0("q=", formatC(as.numeric(gsub("q=", "", stat.test_cytokine_1$p.adj)), format = "e", digits = 1)),
  stat.test_cytokine_1$p.adj
)

bxp_padj <- p +
  stat_pvalue_manual(stat.test_cytokine_1, label = "{q_label}",
                     tip.length = 0.02,
                     inherit.aes = FALSE, hide.ns = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.30)))

ggsave(file.path(FIGURES_DIR, "Figure2D.png"), bxp_padj, dpi = 300, units = "in", width = 6, height = 4)
