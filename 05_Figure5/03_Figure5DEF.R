# ============================================================================
# Purpose:      Figure 5D-F: paired pre vs post vaccination Olink levels of DDX58 (RIG-I), NUB1 and MMP7 in HD/MGUS/SMM, plotted with the plot_cytokine() helper defined in Figure5AB.R.
#               REVISION: q-values are inherited from PANEL_Q, the panel-wide Benjamini-Hochberg family of 156 tests (52 analytes x 3 disease groups) defined in Figure5AB.R, as stated in the Methods.
# Inputs:       Inherited from Figure5AB.R after sourcing it: master_paired data frame, plot_cytokine() and PANEL_Q.
# Outputs:      figures/Figure5D.png, figures/Figure5E.png, figures/Figure5F.png (and matching PDFs).
# Dependencies: R + tidyverse, ggpubr (inherited from Figure5AB.R); sources Figure5AB.R.
# ============================================================================
#!/usr/bin/env Rscript
# Figure 5D-F: DDX58, NUB1, MMP7 paired Olink (pre vs post-vx).
# Reuses data prep and plot_cytokine() from Figure5AB.R.

source("01_Figure5AB.R")  # loads master_paired, plot_cytokine()

# Figure 5D: DDX58 (RIG-I)
p_ddx58 <- plot_cytokine("DDX58", master_paired)
save_figure(p_ddx58, "Figure5D", width = 6, height = 4)

# Figure 5E: NUB1
p_nub1 <- plot_cytokine("NUB1", master_paired)
save_figure(p_nub1, "Figure5E", width = 6, height = 4)

# Figure 5F: MMP7
p_mmp7 <- plot_cytokine("MMP7", master_paired)
save_figure(p_mmp7, "Figure5F", width = 6, height = 4)
