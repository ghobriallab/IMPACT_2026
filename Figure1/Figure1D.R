source("../config.R")
library(tidyverse)
library(lme4)

# Read serial ELISA titer data
dat <- read.csv(file.path(DATA_DIR, "elisa", "elisa_serial_titers_all.csv"))
data_recoded <- dat %>% filter(Days_post2nd >= 14 & Days_post2nd <= 120)
df_swimmer <- data_recoded[, c("Date", "Common_ID", "Disease_Recoded",
                                "Disease_Recoded_granular", "Days_post2nd", "ELISA_Titer")]

# Keep only individuals with serial samples
row_count <- df_swimmer %>%
  group_by(Common_ID) %>%
  summarize(row_count = n()) %>%
  filter(row_count > 1)
df_swimmer <- df_swimmer %>% filter(Common_ID %in% row_count$Common_ID)

# Remove rows where titer increases >0.2 (suggesting re-exposure/boosting)
df_swimmer <- df_swimmer %>%
  group_by(Common_ID) %>%
  mutate(Diff = c(0, diff(ELISA_Titer))) %>%
  mutate(Remove = cumsum(Diff > 0.2)) %>%
  filter(Remove == 0)

# Re-filter for serial samples after removal
row_count <- df_swimmer %>%
  group_by(Common_ID) %>%
  summarize(row_count = n()) %>%
  filter(row_count > 1)
df_swimmer <- df_swimmer %>% filter(Common_ID %in% row_count$Common_ID)

# Order patients by latest sample day
df_swimmer <- df_swimmer %>%
  group_by(Common_ID) %>%
  mutate(max_day = max(na.omit(Days_post2nd))) %>%
  ungroup() %>%
  mutate(Common_ID = fct_reorder(factor(Common_ID), max_day))

# Assign serial sample number per patient
df_swimmer <- df_swimmer %>%
  group_by(Common_ID) %>%
  mutate(ID = factor(row_number())) %>%
  ungroup()

# Subset: untreated SMM and MGUS only; relabel for display
df_swimmer_sub <- df_swimmer %>%
  filter(!Disease_Recoded_granular %in% c("SMM:Treated", "MM", "HD", "IgM-MGUS")) %>%
  mutate(Disease_Recoded_granular = ifelse(Disease_Recoded_granular == "SMM:Untreated", "SMM", Disease_Recoded_granular))

# Calculate per-group counts for facet labels
row_count_sub <- df_swimmer_sub %>%
  group_by(Common_ID, Disease_Recoded_granular) %>%
  summarize(row_count = n(), .groups = "drop")
counts <- row_count_sub %>%
  dplyr::count(Disease_Recoded_granular) %>%
  dplyr::mutate(label = paste(Disease_Recoded_granular, "\n(n =", n, ")"))
labels <- setNames(counts$label, counts$Disease_Recoded_granular)

# Swimmer plot
p1 <- ggplot(df_swimmer_sub, aes(x = Days_post2nd, y = Common_ID, group = Common_ID,
                                  col = Disease_Recoded_granular, shape = factor(ID))) +
  geom_line(linewidth = 0.8, show.legend = TRUE) +
  geom_point(color = "black", fill = "black", size = 2) +
  scale_color_manual(values = c("steelblue", "orange", "tomato2", "chartreuse4"),
                     name = "DiseaseStatus") +
  scale_shape_manual(values = 1:nlevels(df_swimmer$ID), "Serial Sample #") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(-5, 0, 0, 0),
    panel.border = element_rect(fill = NA, color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 10, color = "black")
  ) +
  xlab("Days Post 2nd Dose") + ylab("") +
  guides(fill = "none", color = "none") +
  scale_x_continuous(limits = c(0, NA)) +
  scale_y_discrete(expand = expansion(mult = c(0.03))) +
  facet_wrap(~ Disease_Recoded_granular, scales = "free", ncol = 2,
             labeller = as_labeller(labels))

ggsave(file.path(FIGURES_DIR, "Figure1D.png"), p1, dpi = 300, units = "in",
       height = 4.5, width = 4)
