source("../config.R")
library(tidyverse)

# Read Olink data with tumor burden (M-spike)
bm <- read.csv(file.path(DATA_DIR, "olink", "olink_summary_tumor_burden.csv"),
               check.names = FALSE)

# Filter: Post-vaccination, TNFSF13 (APRIL), SMM patients, QC pass
bm_post <- bm[bm$`Timepoint (red_within 14 days, orange_with pre-vax sample)` == "After_Vax", ]
bm_post <- bm_post[bm_post$Cytokine == "TNFSF13" &
                    bm_post$Diagnosis_when_fully_vaccinated %in% c("SMM") &
                    bm_post$QC_1_pass_0_warning == 1, ]

# Spearman correlation (robust to non-normal distributions)
spearman_test <- cor.test(bm_post$M_spike, bm_post$Level, method = "spearman", exact = FALSE)
rho_value <- signif(spearman_test$estimate, 2)
p_value <- signif(spearman_test$p.value, 2)

p <- ggplot(bm_post) +
  geom_point(aes(M_spike, Level), alpha = 0.75, size = 3, shape = 21,
             color = "black", fill = "steelblue") +
  geom_smooth(aes(M_spike, Level), method = "lm", color = "tomato2",
              fill = "lightblue", level = 0.95) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5)
  ) +
  xlab("M-spike (g/dL)") +
  ylab(expression('Post-Vx ' ~ italic(TNFSF13 ~ (APRIL)) ~ ' level')) +
  ggtitle("SMM patients") +
  annotate("text", x = 3.2, y = 32000,
           label = paste0("Spearman~rho == '", rho_value, "'"),
           parse = TRUE, color = "black", size = 4.5, hjust = 1) +
  annotate("text", x = 3.2, y = 29000,
           label = paste0("p == ", p_value),
           parse = TRUE, color = "black", size = 4.5, hjust = 1)

ggsave(file.path(FIGURES_DIR, "Figure2E.png"), plot = p, dpi = 300, units = "in", width = 3.2, height = 4)
