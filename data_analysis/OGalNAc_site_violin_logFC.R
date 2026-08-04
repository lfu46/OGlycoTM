# OGalNAc_site_violin_logFC.R
# Violin + boxplot of O-GalNAc SITE-level logFC across cell types
# Uses normalized site intensities to compute log2(Tuni/Ctrl) per site

library(tidyverse)
library(ggpubr)
library(rstatix)
library(edgeR)

source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'
figure_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/'

colors_cell <- c("HEK293T" = "#4DBBD5", "HepG2" = "#F39B7F", "Jurkat" = "#00A087")
intensity_cols <- c("Intensity.Tuni_1", "Intensity.Tuni_2", "Intensity.Tuni_3",
                    "Intensity.Ctrl_4", "Intensity.Ctrl_5", "Intensity.Ctrl_6")

# ---- Helper functions (same as O-GlcNAc pipeline) ----
filter_zero_na <- function(data, intensity_cols) {
  data |> filter(if_all(all_of(intensity_cols), ~ . > 0 & !is.na(.)))
}

normalize_intensity <- function(data, intensity_cols) {
  intensity_matrix <- data |> dplyr::select(all_of(intensity_cols)) |> as.matrix()
  col_sums <- colSums(intensity_matrix)
  target_mean <- mean(col_sums)
  sl_factors <- target_mean / col_sums
  intensity_sl <- sweep(intensity_matrix, 2, sl_factors, FUN = "*")
  tmm_factors <- calcNormFactors(intensity_sl)
  intensity_sl_tmm <- sweep(intensity_sl, 2, tmm_factors, FUN = "/")
  intensity_sl_tmm_tb <- as_tibble(intensity_sl_tmm)
  colnames(intensity_sl_tmm_tb) <- paste0(intensity_cols, "_sl_tmm")
  bind_cols(data, intensity_sl_tmm_tb)
}

# ---- Quantify, normalize, compute logFC per site ----
site_logFC_combined <- data.frame()

for (cell in c("HEK293T", "HepG2", "Jurkat")) {
  site <- read_csv(
    paste0(source_file_path, "site/OGalNAc_site_", cell, ".csv"),
    show_col_types = FALSE
  )

  # Aggregate PSMs to site level
  site_quant <- site |>
    group_by(site_index) |>
    summarize(
      Protein.ID = first(Protein.ID),
      Gene = first(Gene),
      across(all_of(intensity_cols), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    )

  # Filter (all channels > 0) and normalize
  site_filtered <- filter_zero_na(site_quant, intensity_cols)
  site_norm <- normalize_intensity(site_filtered, intensity_cols)

  # Compute log2FC per site
  site_fc <- site_norm |>
    mutate(
      mean_Tuni = (Intensity.Tuni_1_sl_tmm + Intensity.Tuni_2_sl_tmm + Intensity.Tuni_3_sl_tmm) / 3,
      mean_Ctrl = (Intensity.Ctrl_4_sl_tmm + Intensity.Ctrl_5_sl_tmm + Intensity.Ctrl_6_sl_tmm) / 3,
      logFC = log2(mean_Tuni / mean_Ctrl),
      CellType = cell
    ) |>
    dplyr::select(site_index, Protein.ID, Gene, logFC, CellType)

  cat(cell, ":", nrow(site_fc), "quantified sites\n")
  site_logFC_combined <- bind_rows(site_logFC_combined, site_fc)
}

site_logFC_combined <- site_logFC_combined |>
  mutate(CellType = factor(CellType, levels = c("HEK293T", "HepG2", "Jurkat")))

# ---- KS tests ----
hek_logfc <- site_logFC_combined |> filter(CellType == "HEK293T") |> pull(logFC)
hepg2_logfc <- site_logFC_combined |> filter(CellType == "HepG2") |> pull(logFC)
jurkat_logfc <- site_logFC_combined |> filter(CellType == "Jurkat") |> pull(logFC)

ks_hh <- ks.test(hek_logfc, hepg2_logfc)
ks_hj <- ks.test(hek_logfc, jurkat_logfc)
ks_gj <- ks.test(hepg2_logfc, jurkat_logfc)

cat("\nKS test: HEK vs HepG2 p =", ks_hh$p.value, "\n")
cat("KS test: HEK vs Jurkat p =", ks_hj$p.value, "\n")
cat("KS test: HepG2 vs Jurkat p =", ks_gj$p.value, "\n")

ks_results <- tibble(
  group1 = c("HEK293T", "HEK293T", "HepG2"),
  group2 = c("HepG2", "Jurkat", "Jurkat"),
  p = c(ks_hh$p.value, ks_hj$p.value, ks_gj$p.value)
) |>
  add_significance("p")

print(ks_results)

# ---- Plot ----
sig_results <- ks_results |> filter(p.signif != "ns")
y_positions <- seq(1.5, by = 0.4, length.out = nrow(sig_results))

Figure6C_site <- site_logFC_combined |>
  ggplot(aes(x = CellType, y = logFC)) +
  geom_violin(aes(fill = CellType), color = "transparent") +
  geom_boxplot(color = "black", outliers = FALSE, width = 0.2, linewidth = 0.3) +
  scale_fill_manual(values = colors_cell) +
  labs(
    x = "",
    y = expression(log[2]*"(Tuni/Ctrl)")
  ) +
  {if (nrow(sig_results) > 0)
    stat_pvalue_manual(
      data = sig_results,
      label = "p.signif",
      tip.length = 0,
      size = 5,
      y.position = y_positions
    )
  } +
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
  filename = paste0(figure_file_path, "Figure6_OGalNAc/Figure6C_OGalNAc_site_logFC_violin.pdf"),
  plot = Figure6C_site,
  width = 1.5, height = 2, units = "in"
)

cat("\nFigure 6C (site-level) saved.\n")
