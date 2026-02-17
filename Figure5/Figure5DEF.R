#!/usr/bin/env Rscript
# Figure 5D-F: DDX58, NUB1, MMP7 paired Olink (pre vs post-vx).
# Reuses data prep and plot_cytokine() from Figure5AB.R.

source("Figure5AB.R")  # loads master_paired, plot_cytokine()

# Figure 5D: DDX58 (RIG-I)
p_ddx58 <- plot_cytokine("DDX58", master_paired)
ggsave(file.path(FIGURES_DIR, "Figure5D.png"), plot = p_ddx58, dpi = 300, units = "in", width = 6, height = 4)

# Figure 5E: NUB1
p_nub1 <- plot_cytokine("NUB1", master_paired)
ggsave(file.path(FIGURES_DIR, "Figure5E.png"), plot = p_nub1, dpi = 300, units = "in", width = 6, height = 4)

# Figure 5F: MMP7
p_mmp7 <- plot_cytokine("MMP7", master_paired)
ggsave(file.path(FIGURES_DIR, "Figure5F.png"), plot = p_mmp7, dpi = 300, units = "in", width = 6, height = 4)
