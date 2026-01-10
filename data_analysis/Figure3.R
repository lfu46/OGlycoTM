# Figure 3: O-GlcNAc Differential Analysis Visualization

library(tidyverse)
library(ggpubr)
library(rstatix)

# Source data paths, colors, and differential analysis results
source('data_source_DE.R')

# Figure 3A ---------------------------------------------------------------
# Violin boxplot comparing O-GlcNAc protein logFC distribution across cell types
# Statistical test: Kolmogorov-Smirnov test (compares distribution shapes)

# Combine logFC results from all cell types
OGlcNAc_protein_logFC_HEK293T <- OGlcNAc_protein_DE_HEK293T %>%
  mutate(CellType = "HEK293T") %>%
  select(Protein.ID, logFC, CellType)

OGlcNAc_protein_logFC_HepG2 <- OGlcNAc_protein_DE_HepG2 %>%
  mutate(CellType = "HepG2") %>%
  select(Protein.ID, logFC, CellType)

OGlcNAc_protein_logFC_Jurkat <- OGlcNAc_protein_DE_Jurkat %>%
  mutate(CellType = "Jurkat") %>%
  select(Protein.ID, logFC, CellType)

OGlcNAc_protein_logFC_combined <- bind_rows(
  OGlcNAc_protein_logFC_HEK293T,
  OGlcNAc_protein_logFC_HepG2,
  OGlcNAc_protein_logFC_Jurkat
)

# Set factor levels for consistent ordering
OGlcNAc_protein_logFC_combined <- OGlcNAc_protein_logFC_combined %>%
  mutate(CellType = factor(CellType, levels = c("HEK293T", "HepG2", "Jurkat")))

# Kolmogorov-Smirnov tests between cell types
ks_HEK293T_HepG2 <- ks.test(
  OGlcNAc_protein_DE_HEK293T$logFC,
  OGlcNAc_protein_DE_HepG2$logFC
)
ks_HEK293T_Jurkat <- ks.test(
  OGlcNAc_protein_DE_HEK293T$logFC,
  OGlcNAc_protein_DE_Jurkat$logFC
)
ks_HepG2_Jurkat <- ks.test(
  OGlcNAc_protein_DE_HepG2$logFC,
  OGlcNAc_protein_DE_Jurkat$logFC
)

# Print KS test results
cat("Figure 3A - KS Test Results:\n")
cat("HEK293T vs HepG2: p =", format(ks_HEK293T_HepG2$p.value, digits = 4), "\n")
cat("HEK293T vs Jurkat: p =", format(ks_HEK293T_Jurkat$p.value, digits = 4), "\n")
cat("HepG2 vs Jurkat: p =", format(ks_HepG2_Jurkat$p.value, digits = 4), "\n")

# Wilcoxon test using rstatix (for comparison)
OGlcNAc_protein_wilcox_test <- OGlcNAc_protein_logFC_combined %>%
  wilcox_test(logFC ~ CellType) %>%
  add_significance("p")

cat("\nFigure 3A - Wilcoxon Test Results:\n")
print(OGlcNAc_protein_wilcox_test)

# Create significance annotation dataframe (using KS test)
OGlcNAc_protein_ks_test <- tribble(
  ~.y., ~group1, ~group2, ~p,
  "logFC", "HEK293T", "HepG2", ks_HEK293T_HepG2$p.value,
  "logFC", "HEK293T", "Jurkat", ks_HEK293T_Jurkat$p.value,
  "logFC", "HepG2", "Jurkat", ks_HepG2_Jurkat$p.value
) %>%
  add_significance("p")

# Create violin boxplot
Figure3A <- OGlcNAc_protein_logFC_combined %>%
  ggplot(aes(x = CellType, y = logFC)) +
  geom_violin(aes(fill = CellType), color = "transparent") +
  geom_boxplot(color = "black", outliers = FALSE, width = 0.2, linewidth = 0.3) +
  scale_fill_manual(values = colors_cell) +
  labs(
    x = "",
    y = expression(log[2]*"(Tuni/Ctrl)")
  ) +
  stat_pvalue_manual(
    data = OGlcNAc_protein_ks_test,
    label = "p.signif",
    tip.length = 0,
    size = 3,
    y.position = c(1.5, 1.9, 1.7)
  ) +
  coord_cartesian(ylim = c(-2, 2.3)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 7),
    axis.text.x = element_text(size = 7, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 7, color = "black"),
    legend.position = "none"
  )

print(Figure3A)

ggsave(
  filename = paste0(figure_file_path, "Figure3/Figure3A.pdf"),
  plot = Figure3A,
  width = 1.5, height = 1.5, units = "in"
)

cat("\nFigure 3A saved to:", figure_file_path, "Figure3/\n")

# Figure 3B ---------------------------------------------------------------
# (Add next figure here)

