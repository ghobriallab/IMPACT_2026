# ============================================================================
# Purpose:      Figure 1E: linear mixed-effects model of antibody-titer waning slope (MGUS vs SMM-untreated), with age and sex as fixed covariates, and the per-individual slope as a random effect.
# Inputs:       data/elisa/elisa_serial_titers_filtered.csv.
# Outputs:      figures/Figure1E.png and the LME fit / LRT statistic for the slope difference.
# Dependencies: R + lme4, lmerTest, dplyr, tidyr, ggplot2; sources ../config.R.
# ============================================================================
source("../config.R")
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(42)

# Configuration ---
INPUT_FILE <- file.path(DATA_DIR, "elisa", "elisa_serial_titers_filtered.csv")
MIN_DAYS <- 14
MAX_DAYS <- 120
MIN_SAMPLES_PER_PATIENT <- 2
DISEASE_GROUPS <- c("MGUS", "SMM:Untreated")

# Load and filter data ---
dat_raw <- read.csv(INPUT_FILE, header = TRUE, stringsAsFactors = FALSE)

dat_filtered <- dat_raw %>%
  filter(Disease_Recoded_granular %in% DISEASE_GROUPS) %>%
  filter(Days_post2nd >= MIN_DAYS & Days_post2nd <= MAX_DAYS) %>%
  mutate(Disease = case_when(
    Disease_Recoded_granular == "MGUS" ~ "MGUS",
    Disease_Recoded_granular == "SMM:Untreated" ~ "SMM",
    TRUE ~ Disease_Recoded_granular
  ))

# Remove samples where titer increases >0.2 (re-exposure/boosting), matching Figure 1D
dat_filtered <- dat_filtered %>%
  arrange(Common_ID, Days_post2nd) %>%
  group_by(Common_ID) %>%
  mutate(Diff = c(0, diff(ELISA_Titer))) %>%
  mutate(Remove = cumsum(Diff > 0.2)) %>%
  filter(Remove == 0) %>%
  select(-Diff, -Remove) %>%
  ungroup()

# age+sex-adjusted waning model. Age at sampling + Sex joined from the deid
# data; restrict to covariate complete-cases BEFORE requiring serial (>= 2) samples so counts/model agree.
dat_filtered <- dat_filtered %>% filter(!is.na(Age) & !is.na(Sex))

# Keep only patients with serial samples (>= 2)
patients_with_serial <- dat_filtered %>%
  group_by(Common_ID) %>%
  summarise(n_samples = n(), .groups = "drop") %>%
  filter(n_samples >= MIN_SAMPLES_PER_PATIENT) %>%
  pull(Common_ID)

dat_serial <- dat_filtered %>% filter(Common_ID %in% patients_with_serial)

samples_by_disease <- dat_serial %>%
  group_by(Disease) %>%
  summarise(
    n_samples = n(),
    n_patients = n_distinct(Common_ID),
    .groups = "drop"
  )

# Mixed-effects models ---
# Log-transform titers for linear modeling
dat_serial <- dat_serial %>% mutate(log_titer = log(ELISA_Titer))

# Full model: disease-specific intercepts AND slopes (random intercept per patient); age+sex covariates
model_full <- lmer(
  log_titer ~ Disease + Days_post2nd:Disease + Age + Sex + (1 | Common_ID),
  data = dat_serial,
  REML = FALSE
)

# Constrained model: disease-specific intercepts, common slope; age+sex covariates
model_constrained <- lmer(
  log_titer ~ Disease + Days_post2nd + Age + Sex + (1 | Common_ID),
  data = dat_serial,
  REML = FALSE
)

# Likelihood ratio test ---
lrt_result <- anova(model_full, model_constrained)
chisq_stat <- lrt_result$Chisq[2]
df_diff <- lrt_result$Df[2] - lrt_result$Df[1]
p_value <- lrt_result$`Pr(>Chisq)`[2]

# Extract disease-specific slopes
coef_summary <- summary(model_full)$coefficients
slopes <- data.frame(
  Disease = c("MGUS", "SMM"),
  Slope = c(
    coef_summary["DiseaseMGUS:Days_post2nd", "Estimate"],
    coef_summary["DiseaseSMM:Days_post2nd", "Estimate"]
  ),
  SE = c(
    coef_summary["DiseaseMGUS:Days_post2nd", "Std. Error"],
    coef_summary["DiseaseSMM:Days_post2nd", "Std. Error"]
  )
)

# Visualization ---
# Predicted regression lines (fixed effects only)
newdata <- expand.grid(
  Days_post2nd = seq(MIN_DAYS, MAX_DAYS, by = 1),
  Disease = c("MGUS", "SMM"),
  Common_ID = "new_patient"
)
# Covariate-adjusted mean trajectory: hold age/sex at cohort mean / modal level for the fixed-effect lines
newdata$Age <- mean(dat_serial$Age, na.rm = TRUE)
newdata$Sex <- names(sort(table(dat_serial$Sex), decreasing = TRUE))[1]
newdata$predicted <- predict(model_full, newdata = newdata, re.form = NA)

# P-values for individual slopes (from lmerTest)
coef_with_pval <- summary(model_full)$coefficients
mgus_pval <- coef_with_pval["DiseaseMGUS:Days_post2nd", "Pr(>|t|)"]
smm_pval <- coef_with_pval["DiseaseSMM:Days_post2nd", "Pr(>|t|)"]

format_pval <- function(p) {
  # show the actual p value (scientific notation when small)
  if (p < 0.001) return(sprintf("p=%.1e", p))
  else return(sprintf("p=%.3f", p))
}

# Annotation text
annot_mgus <- sprintf("MGUS: \u03B2=%.4f/day, %s", slopes$Slope[1], format_pval(mgus_pval))
annot_smm <- sprintf("SMM: \u03B2=%.4f/day, %s", slopes$Slope[2], format_pval(smm_pval))
annot_lrt <- sprintf("LRT %s", format_pval(p_value))

y_min <- min(dat_serial$log_titer)
y_range <- max(dat_serial$log_titer) - y_min
annot_y_base <- y_min + 0.08 * y_range

p <- ggplot() +
  geom_line(data = dat_serial,
            aes(x = Days_post2nd, y = log_titer, group = Common_ID, color = Disease),
            alpha = 0.4, linewidth = 0.4) +
  geom_point(data = dat_serial,
             aes(x = Days_post2nd, y = log_titer, color = Disease),
             alpha = 0.5, size = 1) +
  geom_line(data = newdata,
            aes(x = Days_post2nd, y = predicted, color = Disease),
            linewidth = 1.2) +
  scale_color_manual(
    values = c("MGUS" = "steelblue", "SMM" = "tomato"),
    labels = c(
      "MGUS" = paste0("MGUS (n=", samples_by_disease$n_patients[samples_by_disease$Disease == "MGUS"], ")"),
      "SMM" = paste0("SMM (n=", samples_by_disease$n_patients[samples_by_disease$Disease == "SMM"], ")")
    )
  ) +
  annotate("text", x = 5, y = annot_y_base + 0.40,
           label = annot_mgus, hjust = 0, size = 2.5, color = "steelblue", fontface = "bold") +
  annotate("text", x = 5, y = annot_y_base + 0.18,
           label = annot_smm, hjust = 0, size = 2.5, color = "tomato", fontface = "bold") +
  annotate("text", x = 5, y = annot_y_base - 0.04,
           label = annot_lrt, hjust = 0, size = 2.5, color = "black", fontface = "bold") +
  xlab("Days post 2nd dose") +
  ylab("log(ELISA Titer)") +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    legend.title = element_blank(),
    legend.position = c(0.02, 0.38),
    legend.justification = c(0, 0),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.background = element_blank(),
    legend.key = element_blank(),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7, color = "black"),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    plot.margin = margin(5, 5, 5, 5)
  )

ggsave(file.path(FIGURES_DIR, "Figure1E.png"), p, width = 3.0, height = 2.8, dpi = 300)
