source("../config.R")
library(tidyverse)
library(rstatix)
library(ggpubr)

# Read Olink proteomics data
data <- read_csv(file.path(DATA_DIR, "olink", "olink_cytokines.csv"))

# Keep only paired samples (pre + post vaccination)
plot_df_index <- data %>%
  group_by(Patient_ID, Cytokine) %>%
  filter(any(Timepoint == "Pre-Vx"), any(Timepoint == "Post-Vx"))

plot_df_index$Timepoint <- ifelse(plot_df_index$Timepoint == "Pre-Vx", "1", "2")
plot_df_index$Timepoint <- factor(plot_df_index$Timepoint, levels = c("1", "2"), ordered = TRUE)
plot_df_index$dummy_timepoint <- as.numeric(plot_df_index$Timepoint)

# Add jitter for paired lines
b <- as.list(runif(nrow(plot_df_index), -0.15, 0.15))
b_df <- do.call("rbind", b)
colnames(b_df) <- "add_noise"
plot_df_index <- cbind(plot_df_index, plot_df_index$dummy_timepoint + b_df)
colnames(plot_df_index)[ncol(plot_df_index)] <- "add_noise"

# Ensure paired filtering
master_paired <- plot_df_index %>%
  group_by(Patient_ID) %>%
  filter(any(Timepoint == "1"), any(Timepoint == "2"))

# Build facet labels with sample counts
counts <- as.data.frame(table(master_paired$Timepoint, master_paired$Diagnosis, master_paired$Cytokine))
counts <- counts %>% distinct(Var2, Var3, .keep_all = TRUE)
counts$matchname <- paste0(counts$Var2, "_", counts$Var3)
master_paired$matchname <- paste0(master_paired$Diagnosis, "_", master_paired$Cytokine)
counts$n_ind <- paste0(counts$Var2, "\n (n=", counts$Freq, ")")
master_paired$n_ind <- counts$n_ind[match(master_paired$matchname, counts$matchname)]

# Display names for y-axis labels
cyto_labels <- c("IL1B" = "IL-1\u03B2", "IL18" = "IL-18",
                 "NUB1" = "NUB1", "DDX58" = "DDX58", "MMP7" = "MMP7")

# --- Helper function for individual cytokine plots ---
plot_cytokine <- function(cyto_name, paired_data) {
  plot_paired_df <- paired_data %>% filter(Cytokine %in% cyto_name)
  plot_paired_df$VaccineTimepoint <- plot_paired_df$Timepoint
  y_label <- ifelse(cyto_name %in% names(cyto_labels), cyto_labels[cyto_name], cyto_name)

  p <- ggplot(plot_paired_df, aes(x = VaccineTimepoint, y = Level, fill = VaccineTimepoint)) +
    geom_violin(alpha = 0.5, position = position_dodge(width = .75), size = 1, show.legend = FALSE, color = NA) +
    geom_boxplot(alpha = 0.5, show.legend = FALSE) +
    geom_line(aes(add_noise, y = Level, group = Patient_ID),
              color = "black", size = 0.3, alpha = 0.5, show.legend = FALSE) +
    geom_pointrange(aes(add_noise, ymin = Level, ymax = Level),
                    shape = 21, size = 0.5, alpha = 1, show.legend = FALSE, linetype = "dotted") +
    scale_fill_manual(values = c("steelblue", "tomato2")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(fill = NA, color = "black"),
          axis.text = element_text(size = 16, color = "black"),
          axis.title = element_text(size = 18, color = "black"),
          strip.text = element_text(size = 18, face = "plain"),
          strip.background = element_blank()) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    scale_x_discrete(breaks = c(1, 2), labels = c("Pre-Vx", "Post-Vx")) +
    facet_wrap(. ~ n_ind) +
    ylab(paste0(y_label, " (NPX)")) + xlab("")

  # Paired Wilcoxon test per disease group, BH-corrected
  stat.test <- plot_paired_df %>%
    group_by(n_ind) %>%
    wilcox_test(Level ~ VaccineTimepoint, paired = TRUE) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance()

  # Paired rank-biserial effect sizes per disease group
  eff <- plot_paired_df %>%
    group_by(n_ind) %>%
    wilcox_effsize(Level ~ VaccineTimepoint, paired = TRUE)

  stat.test <- stat.test %>%
    add_xy_position(x = "VaccineTimepoint") %>%
    left_join(eff %>% select(n_ind, effsize), by = "n_ind") %>%
    mutate(
      q_num = signif(p.adj, 2),
      q_fmt = paste0("q=", q_num),
      bracket_label = ifelse(q_num < 0.1,
                             paste0(q_fmt, ", r=", round(abs(effsize), 2)),
                             q_fmt))

  final <- p +
    stat_pvalue_manual(stat.test, label = "bracket_label",
                       tip.length = 0.02, bracket.nudge.y = 1, size = 5,
                       inherit.aes = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

  return(final)
}

# Figure 5A: IL-1B
p_il1b <- plot_cytokine("IL1B", master_paired)
ggsave(file.path(FIGURES_DIR, "Figure5A.png"), plot = p_il1b, dpi = 300, units = "in", width = 6, height = 4)

# Figure 5B: IL-18
p_il18 <- plot_cytokine("IL18", master_paired)
ggsave(file.path(FIGURES_DIR, "Figure5B.png"), plot = p_il18, dpi = 300, units = "in", width = 6, height = 4)
