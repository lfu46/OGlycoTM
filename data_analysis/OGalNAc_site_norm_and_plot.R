# OGalNAc_site_norm_and_plot.R
# 1. Quantify O-GalNAc sites (aggregate PSMs per site)
# 2. Normalize (SL + TMM)
# 3. Plot ranking with normalized fold change

library(tidyverse)
library(edgeR)

source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'
figure_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/'

colors_cell <- c("HEK293T" = "#4DBBD5", "HepG2" = "#F39B7F", "Jurkat" = "#00A087")
intensity_cols <- c("Intensity.Tuni_1", "Intensity.Tuni_2", "Intensity.Tuni_3",
                    "Intensity.Ctrl_4", "Intensity.Ctrl_5", "Intensity.Ctrl_6")

dir.create(paste0(figure_file_path, "Figure5_OGalNAc"), showWarnings = FALSE)

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
  intensity_sl_tb <- as_tibble(intensity_sl)
  colnames(intensity_sl_tb) <- paste0(intensity_cols, "_sl")
  tmm_factors <- calcNormFactors(intensity_sl)
  intensity_sl_tmm <- sweep(intensity_sl, 2, tmm_factors, FUN = "/")
  intensity_sl_tmm_tb <- as_tibble(intensity_sl_tmm)
  colnames(intensity_sl_tmm_tb) <- paste0(intensity_cols, "_sl_tmm")
  bind_cols(data, intensity_sl_tb, intensity_sl_tmm_tb)
}

# ---- Process each cell type ----

for (cell in c("HEK293T", "HepG2", "Jurkat")) {
  cat("\n", strrep("=", 40), "\n", cell, "\n", strrep("=", 40), "\n")

  site <- read_csv(
    paste0(source_file_path, "site/OGalNAc_site_", cell, ".csv"),
    show_col_types = FALSE
  )

  # Quantify: sum PSM intensities per site
  site_quant <- site |>
    group_by(site_index) |>
    summarize(
      Protein.ID = first(Protein.ID),
      Gene = first(Gene),
      Protein.Description = first(Protein.Description),
      modified_residue = first(modified_residue),
      site_number = first(site_number),
      across(all_of(intensity_cols), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  cat("  Quantified sites:", nrow(site_quant), "\n")

  write_csv(site_quant, paste0(source_file_path, "quantification/OGalNAc_site_quant_", cell, ".csv"))

  # Filter and normalize
  site_filtered <- filter_zero_na(site_quant, intensity_cols)
  cat("  After filtering:", nrow(site_filtered), "\n")

  site_norm <- normalize_intensity(site_filtered, intensity_cols)

  write_csv(site_norm, paste0(source_file_path, "normalization/OGalNAc_site_norm_", cell, ".csv"))

  # Calculate log2FC from normalized data
  site_plot <- site_norm |>
    mutate(
      mean_Tuni = (Intensity.Tuni_1_sl_tmm + Intensity.Tuni_2_sl_tmm + Intensity.Tuni_3_sl_tmm) / 3,
      mean_Ctrl = (Intensity.Ctrl_4_sl_tmm + Intensity.Ctrl_5_sl_tmm + Intensity.Ctrl_6_sl_tmm) / 3,
      log2FC = log2(mean_Tuni / mean_Ctrl)
    ) |>
    arrange(log2FC) |>
    mutate(rank = row_number())

  cat("  log2FC range:", round(min(site_plot$log2FC), 2), "to", round(max(site_plot$log2FC), 2), "\n")
  cat("  Up (>0.5):", sum(site_plot$log2FC > 0.5), ", Down (<-0.5):", sum(site_plot$log2FC < -0.5), "\n")

  # Plot
  p <- site_plot |>
    ggplot(aes(x = rank, y = log2FC)) +
    geom_point(size = 0.8, color = colors_cell[cell], alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    labs(
      x = "O-GalNAc site rank",
      y = expression(log[2](Tuni/Ctrl))
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 9, color = "black"),
      axis.text = element_text(size = 9, color = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )

  ggsave(
    paste0(figure_file_path, "Figure5_OGalNAc/OGalNAc_site_ranking_", cell, ".pdf"),
    p, width = 2, height = 1.5, units = "in"
  )
  cat("  Plot saved\n")
}

cat("\nAll done.\n")
