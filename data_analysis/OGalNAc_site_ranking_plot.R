# OGalNAc_site_ranking_plot.R
# Site-level ranking dot plot for O-GalNAc glycopeptides
# Each dot = one localized site, ranked by log2(Tuni/Ctrl)

library(tidyverse)

source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'
figure_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/'

colors_cell <- c("HEK293T" = "#4DBBD5", "HepG2" = "#F39B7F", "Jurkat" = "#00A087")

dir.create(paste0(figure_file_path, "Figure5_OGalNAc"), showWarnings = FALSE)

for (cell in c("HEK293T", "HepG2", "Jurkat")) {
  site <- read_csv(
    paste0(source_file_path, "site/OGalNAc_site_", cell, ".csv"),
    show_col_types = FALSE
  )

  # Aggregate to unique site (sum intensities per site_index)
  site_quant <- site |>
    group_by(site_index) |>
    summarize(
      Gene = first(Gene),
      Intensity.Tuni_1 = sum(Intensity.Tuni_1, na.rm = TRUE),
      Intensity.Tuni_2 = sum(Intensity.Tuni_2, na.rm = TRUE),
      Intensity.Tuni_3 = sum(Intensity.Tuni_3, na.rm = TRUE),
      Intensity.Ctrl_4 = sum(Intensity.Ctrl_4, na.rm = TRUE),
      Intensity.Ctrl_5 = sum(Intensity.Ctrl_5, na.rm = TRUE),
      Intensity.Ctrl_6 = sum(Intensity.Ctrl_6, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      mean_Tuni = (Intensity.Tuni_1 + Intensity.Tuni_2 + Intensity.Tuni_3) / 3,
      mean_Ctrl = (Intensity.Ctrl_4 + Intensity.Ctrl_5 + Intensity.Ctrl_6) / 3
    ) |>
    filter(mean_Ctrl > 0 & mean_Tuni > 0) |>
    mutate(log2FC = log2(mean_Tuni / mean_Ctrl)) |>
    arrange(log2FC) |>
    mutate(rank = row_number())

  n_sites <- nrow(site_quant)
  cat(cell, ":", n_sites, "quantified sites\n")

  p <- site_quant |>
    ggplot(aes(x = rank, y = log2FC)) +
    geom_point(size = 1.2, color = colors_cell[cell], alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    labs(
      x = "Glycopeptide rank",
      y = expression(log[2](Tuni/Ctrl))
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9, color = "black"),
      plot.margin = margin(5, 10, 5, 5)
    )

  ggsave(
    paste0(figure_file_path, "Figure5_OGalNAc/OGalNAc_site_ranking_", cell, ".pdf"),
    p, width = 2.5, height = 2, units = "in"
  )
  cat("  Saved:", cell, "\n")
}

cat("\nAll site ranking plots saved.\n")
