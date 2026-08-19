# ============================================================================
# Purpose:      Figure 2D: panel-wide Olink screen. Every protein on the panel is tested for a
#               disease-vs-HD difference at both timepoints, and the resulting q-values are plotted
#               against each other, so the reader can see that APRIL (TNFSF13) is not a pre-selected
#               candidate but the top hit of an unbiased screen. All three pairwise contrasts are
#               shown, SMM vs HD, MGUS vs HD and SMM vs MGUS, so the panel displays the same 156
#               tests the q-values were corrected across. Comparisons are age- and sex-adjusted rank-based ANCOVA, Benjamini-
#               Hochberg corrected across all 156 tests within a timepoint (52 proteins x 3 pairwise
#               contrasts), the same family used throughout the Olink analysis.
# Inputs:       data/olink/olink_paired_prepost.csv, data/elisa/elisa_cohort_demographics.csv.
# Outputs:      figures/Figure2D.png, .pdf and .svg; source data
#               tables/Figure2D_Olink_screen_PreVx.csv and tables/Figure2D_Olink_screen_PostVx.csv.
# Dependencies: R + tidyverse, rstatix, ggrepel; sources ../config.R for FONT and save_figure().
# ============================================================================
source("../config.R")
library(tidyverse)
library(rstatix)
library(ggrepel)

# ---------------------------------------------------------------- data, as in the other Fig 2 panels
olink <- read.csv(file.path(DATA_DIR, "olink", "olink_paired_prepost.csv"), check.names = FALSE)
olink$Diagnosis <- olink$Diagnosis_when_fully_vaccinated
olink_meta <- olink[!duplicated(olink$Patient_ID), ]

olink <- olink[!is.na(olink$Level), ]
olink <- olink[olink$Diagnosis != "IgM-MGUS", ]

# The single MGUS participant with prior systemic therapy is dropped so the comparison stays
# treatment-naive, matching the other Olink panels.
.OLINK_TREATED_DROP <- c("P49")
olink <- olink[!olink$Patient_ID %in% .OLINK_TREATED_DROP, ]

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

data <- olink_2 %>% filter(QC_1_pass_0_warning == "1")
data <- data[data$Diagnosis != "IgM-MGUS", ]

# Proteins with more than one measurement in every group, which yields the 52-analyte panel
check_data <- as.data.frame(table(data$Cytokine, data$Diagnosis, data$Timepoint))
colnames(check_data) <- c("Cytokine", "Diagnosis", "Timepoint", "Freq")
check_data_wide <- reshape(check_data, idvar = c("Cytokine", "Timepoint"),
                           timevar = "Diagnosis", direction = "wide")
cytokines_with_enough_data <- check_data_wide %>%
  group_by(Timepoint) %>%
  filter(Freq.Healthy > 1 & Freq.MGUS > 1 & Freq.SMM > 1) %>%
  pull(Cytokine)
cytokines_with_enough_data <- unique(as.character(cytokines_with_enough_data))

data$Timepoint <- ifelse(data$Timepoint == "After_Vax", "Post-Vx", "Pre-Vx")
data <- data %>% filter(Cytokine %in% cytokines_with_enough_data)

# ---------------------------------------------------------------- covariates
# Age and sex are taken from the de-identified ELISA cohort table, NOT from the Olink file. The Age
# column shipped with the Olink data carries a single fill value for every healthy donor, which
# would make age a proxy for group membership and absorb the disease effect. The assertion below
# fails loudly if a constant-age source is ever reintroduced.
demo <- read.csv(file.path(DATA_DIR, "elisa", "elisa_cohort_demographics.csv")) %>%
  filter(!is.na(Age_at_second_dose), !is.na(Sex)) %>%
  distinct(Deidentified_Patient_ID, .keep_all = TRUE) %>%
  transmute(Patient_ID = as.character(Deidentified_Patient_ID),
            Age = Age_at_second_dose, Sex = Sex)
data <- data %>% left_join(demo, by = "Patient_ID")
stopifnot(!any(is.na(data$Age)), !any(is.na(data$Sex)))
stopifnot(dplyr::n_distinct(data$Age[data$Diagnosis == "Healthy"]) > 1)

# ---------------------------------------------------------------- the screen
# Age- and sex-adjusted rank-based ANCOVA, the covariate-adjusted analogue of the Wilcoxon rank-sum
# test: the rank-transformed outcome is regressed on disease group, age and sex, and the group
# coefficient is tested. The unadjusted Wilcoxon q is carried alongside in the source data so either
# can be checked.
PAIRS <- list(c("Healthy", "MGUS"), c("Healthy", "SMM"), c("MGUS", "SMM"))

screen_one <- function(tp, cyto, pr) {
  s <- data %>% filter(Timepoint == tp, Cytokine == cyto, Diagnosis %in% pr)
  if (dplyr::n_distinct(s$Diagnosis) < 2) return(NULL)
  s$rk <- rank(s$Level)
  cf <- coef(summary(lm(rk ~ Diagnosis + Age + factor(Sex), data = s)))
  # Same call the Figure 2E panel makes, so the unadjusted column reproduces the published
  # q-values exactly. rstatix uses the exact test at these sample sizes; forcing the normal
  # approximation here would silently disagree with the rest of the paper.
  wt <- rstatix::wilcox_test(s, Level ~ Diagnosis)
  m  <- tapply(s$Level, s$Diagnosis, median)
  tibble(Timepoint = tp, Protein = cyto, Contrast = paste(pr, collapse = " vs "),
         n_group1 = sum(s$Diagnosis == pr[1]), n_group2 = sum(s$Diagnosis == pr[2]),
         median_group1 = unname(m[pr[1]]), median_group2 = unname(m[pr[2]]),
         direction = ifelse(m[pr[2]] < m[pr[1]], paste("lower in", pr[2]),
                            paste("higher in", pr[2])),
         p_adjusted = cf[2, 4], p_unadjusted = wt$p[1])
}

res <- map_dfr(c("Pre-Vx", "Post-Vx"), function(tp)
  map_dfr(cytokines_with_enough_data, function(cy)
    map_dfr(PAIRS, function(pr) screen_one(tp, cy, pr)))) %>%
  group_by(Timepoint) %>%
  mutate(q_adjusted   = p.adjust(p_adjusted,   method = "BH"),
         q_unadjusted = p.adjust(p_unadjusted, method = "BH")) %>%
  ungroup()

stopifnot(length(cytokines_with_enough_data) == 52)
stopifnot(all(table(res$Timepoint) == 156))
message(sprintf("panel-wide screen: %d proteins x %d contrasts = %d tests per timepoint",
                length(cytokines_with_enough_data), length(PAIRS), 156))

# source data, one file per timepoint, ranked by the adjusted q
for (tp in c("Pre-Vx", "Post-Vx")) {
  out <- res %>% filter(Timepoint == tp) %>% arrange(q_adjusted, p_adjusted) %>%
    mutate(rank = row_number(),
           across(c(median_group1, median_group2), ~ round(.x, 3)),
           across(c(p_adjusted, q_adjusted, p_unadjusted, q_unadjusted), ~ signif(.x, 4))) %>%
    select(rank, Protein, Contrast, n_group1, n_group2, median_group1, median_group2,
           direction, p_adjusted, q_adjusted, p_unadjusted, q_unadjusted)
  dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)
  write.csv(out, file.path(TABLES_DIR, sprintf("Figure2D_Olink_screen_%s.csv",
                                               sub("-", "", tp))), row.names = FALSE)
}

# ---------------------------------------------------------------- panel
Q_THRESH <- 0.1
# All three contrasts of the correction family are drawn, so the panel shows the same 156 tests
# the q-values were corrected across rather than a subset of them.
LAB <- c("Healthy vs SMM" = "SMM vs HD", "Healthy vs MGUS" = "MGUS vs HD",
         "MGUS vs SMM" = "SMM vs MGUS")

plot_df <- res %>%
  filter(Contrast %in% names(LAB)) %>%
  select(Protein, Contrast, Timepoint, q_adjusted) %>%
  pivot_wider(names_from = Timepoint, values_from = q_adjusted) %>%
  mutate(x = -log10(`Pre-Vx`), y = -log10(`Post-Vx`),
         Comparison = factor(LAB[Contrast],
                             levels = c("SMM vs HD", "MGUS vs HD", "SMM vs MGUS")),
         sig_both = `Pre-Vx` < Q_THRESH & `Post-Vx` < Q_THRESH)

# Benjamini-Hochberg is a step-up procedure, so many proteins inherit exactly the same q and their
# points would sit on top of one another: 15 of them share a single coordinate without this. The
# offset is applied once, from a fixed seed, and stored as a column so the point and its label use
# the same position. It is well below the distance to the significance lines, so no point can be
# nudged across a threshold.
JITTER <- 0.018
set.seed(2026)
plot_df <- plot_df %>%
  mutate(xj = x + runif(n(), -JITTER, JITTER),
         yj = y + runif(n(), -JITTER, JITTER))
message(sprintf("plotted points: %d (%d proteins x %d contrasts)",
                nrow(plot_df), n_distinct(plot_df$Protein), length(LAB)))

# Only APRIL is labelled; every other protein is identifiable from the source-data CSVs.
lab_df <- plot_df %>% filter(Protein == "TNFSF13", Comparison == "SMM vs HD")
# The two axes span very different ranges (the strongest pre-Vx hit is far weaker than
# the strongest post-Vx hit), so a single square limit left roughly 40% of the panel empty. Each
# axis now gets its own data-driven limit, which fills the wider canvas.
xlim_hi <- max(plot_df$xj) * 1.13
ylim_hi <- max(plot_df$yj) * 1.10

p <- ggplot(plot_df, aes(xj, yj, fill = Comparison)) +
  geom_hline(yintercept = -log10(Q_THRESH), linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = -log10(Q_THRESH), linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  geom_point(shape = 21, size = 3.3, colour = "black", stroke = 0.4, alpha = 0.85) +
  geom_text_repel(data = lab_df, aes(label = "TNFSF13 (APRIL)"),
                  size = 4.6, fontface = "italic", colour = "black", family = FONT,
                  min.segment.length = 0, box.padding = 0.9, seed = 2026,
                  segment.colour = "grey30", segment.size = 0.3, show.legend = FALSE) +
  annotate("text", x = xlim_hi * 0.98, y = -log10(Q_THRESH), label = "q = 0.1",
           hjust = 1, vjust = -0.5, size = 3.8, colour = "grey35", family = FONT) +
  # MGUS takes the same mustard/orange used for MGUS in the APRIL panel (03_Figure2E.R);
  # SMM keeps the project red.
  scale_fill_manual(values = c("SMM vs HD" = unname(COLORS["SMM"]),
                               "MGUS vs HD" = "orange",
                               "SMM vs MGUS" = "#3C5488")) +
  coord_cartesian(xlim = c(0, xlim_hi), ylim = c(0, ylim_hi)) +
  labs(title = "Olink screen, 52 proteins",
       subtitle = paste("Age- and sex-adjusted rank-based ANCOVA\n",
                        "All 156 tests per timepoint, BH-corrected together; points jittered",
                        sep = ""),
       # Plain Unicode rather than plotmath: the plotmath unary minus collides with the "l" of log
       # in the default device font.
       x = "Pre-Vx  \u2212log\u2081\u2080(q)",
       y = "Post-Vx  \u2212log\u2081\u2080(q)", fill = NULL) +
  theme_bw(base_family = FONT) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
        panel.border = element_rect(fill = NA, colour = "black"),
        text = element_text(family = FONT),
        axis.text = element_text(size = 13, colour = "black", family = FONT),
        axis.title = element_text(size = 15, family = FONT),
        # bottom-right is the empty corner: the null cloud sits bottom-left and the only
        # significant point sits top-right
        legend.position = c(0.985, 0.015), legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
        legend.key.size = unit(0.85, "lines"),
        legend.key.height = unit(1.05, "lines"),
        legend.text = element_text(size = 13, family = FONT),
        plot.title = element_text(size = 15, family = FONT, hjust = 0,
                                  margin = margin(b = 1)),
        plot.subtitle = element_text(size = 11, family = FONT, colour = "grey25",
                                     hjust = 0, lineheight = 1.1,
                                     margin = margin(t = 0, b = 4)),
        plot.margin = margin(5, 7, 3, 3))

save_figure(p, "Figure2D", width = 7.0, height = 3.85)

# console summary
cat("\nproteins clearing q<", Q_THRESH, " at BOTH timepoints:\n", sep = "")
print(as.data.frame(plot_df %>% filter(sig_both) %>%
  select(Protein, Comparison, `Pre-Vx`, `Post-Vx`) %>% arrange(`Post-Vx`)))
cat("\nany-timepoint hits per comparison:\n")
print(plot_df %>% group_by(Comparison) %>%
  summarise(n_proteins = n(),
            sig_pre = sum(`Pre-Vx` < Q_THRESH), sig_post = sum(`Post-Vx` < Q_THRESH),
            sig_both = sum(sig_both), .groups = "drop"))
