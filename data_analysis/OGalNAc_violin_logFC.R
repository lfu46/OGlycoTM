# OGalNAc_violin_logFC.R
# Violin + boxplot of O-GalNAc protein logFC across cell types
# Matches Figure 3A style exactly

library(tidyverse)
library(ggpubr)
library(rstatix)

source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'
figure_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/'

colors_cell <- c("HEK293T" = "#4DBBD5", "HepG2" = "#F39B7F", "Jurkat" = "#00A087")

# Load O-GalNAc DE results
OGalNAc_logFC_combined <- bind_rows(
  read_csv(paste0(source_file_path, "differential_analysis/OGalNAc_protein_DE_HEK293T.csv"), show_col_types = FALSE) |>
    mutate(CellType = "HEK293T") |> dplyr::select(Protein.ID, logFC, CellType),
  read_csv(paste0(source_file_path, "differential_analysis/OGalNAc_protein_DE_HepG2.csv"), show_col_types = FALSE) |>
    mutate(CellType = "HepG2") |> dplyr::select(Protein.ID, logFC, CellType),
  read_csv(paste0(source_file_path, "differential_analysis/OGalNAc_protein_DE_Jurkat.csv"), show_col_types = FALSE) |>
    mutate(CellType = "Jurkat") |> dplyr::select(Protein.ID, logFC, CellType)
) |>
  mutate(CellType = factor(CellType, levels = c("HEK293T", "HepG2", "Jurkat")))

# KS tests between cell types (same as Figure 3A)
hek_logfc <- OGalNAc_logFC_combined |> filter(CellType == "HEK293T") |> pull(logFC)
hepg2_logfc <- OGalNAc_logFC_combined |> filter(CellType == "HepG2") |> pull(logFC)
jurkat_logfc <- OGalNAc_logFC_combined |> filter(CellType == "Jurkat") |> pull(logFC)

ks_hh <- ks.test(hek_logfc, hepg2_logfc)
ks_hj <- ks.test(hek_logfc, jurkat_logfc)
ks_gj <- ks.test(hepg2_logfc, jurkat_logfc)

cat("KS test: HEK vs HepG2 p =", ks_hh$p.value, "\n")
cat("KS test: HEK vs Jurkat p =", ks_hj$p.value, "\n")
cat("KS test: HepG2 vs Jurkat p =", ks_gj$p.value, "\n")

ks_results <- tibble(
  group1 = c("HEK293T", "HEK293T", "HepG2"),
  group2 = c("HepG2", "Jurkat", "Jurkat"),
  p = c(ks_hh$p.value, ks_hj$p.value, ks_gj$p.value)
) |>
  add_significance("p")

# Plot
Figure6E <- OGalNAc_logFC_combined |>
  ggplot(aes(x = CellType, y = logFC)) +
  geom_violin(aes(fill = CellType), color = "transparent") +
  geom_boxplot(color = "black", outliers = FALSE, width = 0.2, linewidth = 0.3) +
  scale_fill_manual(values = colors_cell) +
  labs(
    x = "",
    y = expression(log[2]*"(Tuni/Ctrl)")
  ) +
  stat_pvalue_manual(
    data = ks_results |> filter(p.signif != "ns"),
    label = "p.signif",
    tip.length = 0,
    size = 5,
    y.position = c(1.5)
  ) +
  coord_cartesian(ylim = c(-2, 2.3)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 9, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.position = "none"
  )

ggsave(
  filename = paste0(figure_file_path, "Figure6_OGalNAc/Figure6C_OGalNAc_logFC_violin.pdf"),
  plot = Figure6E,
  width = 1.5, height = 2, units = "in"
)

cat("Figure 6E saved.\n")
