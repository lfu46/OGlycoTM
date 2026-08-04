# =============================================================================
# Figure 4C Candidate Plots for Slide Review
# =============================================================================
# Generate one dot plot per candidate protein in the exact Figure 4C style.
# Candidates are shown in every cell type where they were quantified (no
# truncation). A single fixed y-axis range is used across all candidates so
# the slides are visually comparable side-by-side.
#
# Candidates (gene -> UniProt):
#   CCAR1 -> Q8IX12  (3 cells)
#   YIPF3 -> Q9GZM5  (HepG2 + Jurkat only)
#   CIC   -> Q96RK0  (3 cells)
#   YIF1B -> Q5BJH7  (3 cells)

library(tidyverse)

source('data_source.R')

out_dir <- paste0(figure_file_path, 'Figure4/Figure4C_candidate_slides/')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load normalized O-GlcNAc protein data --------------------------------
OGlcNAc_protein_norm_HEK293T <- read_csv(
  paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_HEK293T.csv'),
  show_col_types = FALSE
)
OGlcNAc_protein_norm_HepG2 <- read_csv(
  paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_HepG2.csv'),
  show_col_types = FALSE
)
OGlcNAc_protein_norm_Jurkat <- read_csv(
  paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_Jurkat.csv'),
  show_col_types = FALSE
)

# ---- Candidate table ------------------------------------------------------
candidates <- tibble::tribble(
  ~Gene,    ~Protein_ID,
  "CCAR1",  "Q8IX12",
  "YIPF3",  "Q9GZM5",
  "CIC",    "Q96RK0",
  "YIF1B",  "Q5BJH7"
)

# ---- Per-cell fold-change helper (same logic as existing Figure 4C) -------
calculate_fc <- function(norm_data, cell_name, protein_id) {
  norm_data |>
    filter(Protein.ID == protein_id) |>
    dplyr::select(Protein.ID, Gene,
                  Intensity.Tuni_1_sl_tmm, Intensity.Tuni_2_sl_tmm, Intensity.Tuni_3_sl_tmm,
                  Intensity.Ctrl_4_sl_tmm, Intensity.Ctrl_5_sl_tmm, Intensity.Ctrl_6_sl_tmm) |>
    rowwise() |>
    mutate(
      mean_Ctrl   = mean(c(Intensity.Ctrl_4_sl_tmm, Intensity.Ctrl_5_sl_tmm, Intensity.Ctrl_6_sl_tmm), na.rm = TRUE),
      log2FC_rep1 = log2(Intensity.Tuni_1_sl_tmm / mean_Ctrl),
      log2FC_rep2 = log2(Intensity.Tuni_2_sl_tmm / mean_Ctrl),
      log2FC_rep3 = log2(Intensity.Tuni_3_sl_tmm / mean_Ctrl)
    ) |>
    ungroup() |>
    dplyr::select(Protein.ID, Gene, log2FC_rep1, log2FC_rep2, log2FC_rep3) |>
    pivot_longer(cols = starts_with("log2FC"),
                 names_to = "replicate",
                 values_to = "log2FC") |>
    mutate(cell = cell_name)
}

# ---- Build combined dataset for all candidates ----------------------------
all_fc <- purrr::map_dfr(seq_len(nrow(candidates)), function(i) {
  pid  <- candidates$Protein_ID[i]
  gene <- candidates$Gene[i]
  bind_rows(
    calculate_fc(OGlcNAc_protein_norm_HEK293T, "HEK293T", pid),
    calculate_fc(OGlcNAc_protein_norm_HepG2,   "HepG2",   pid),
    calculate_fc(OGlcNAc_protein_norm_Jurkat,  "Jurkat",  pid)
  ) |>
    mutate(Gene_target = gene)
}) |>
  filter(!is.na(log2FC)) |>
  mutate(cell = factor(cell, levels = c("HEK293T", "HepG2", "Jurkat")))

cat("\nCombined candidate data:\n")
print(all_fc |> group_by(Gene_target, cell) |> summarise(n = n(), .groups = "drop"))

# ---- Fixed y-axis across candidates ---------------------------------------
y_min <- floor(min(all_fc$log2FC, na.rm = TRUE) * 2) / 2
y_max <- ceiling(max(all_fc$log2FC, na.rm = TRUE) * 2) / 2
# pad slightly so points aren't clipped
y_lims <- c(y_min - 0.1, y_max + 0.1)
cat(sprintf("\nFixed y-axis limits: [%.2f, %.2f]\n", y_lims[1], y_lims[2]))

# ---- Plot helper (matches Figure 4C style) --------------------------------
plot_candidate <- function(gene_name, protein_id) {
  df <- all_fc |>
    filter(Gene_target == gene_name) |>
    mutate(cell = droplevels(cell))

  p <- df |>
    ggplot(aes(x = cell, y = log2FC, color = cell)) +
    geom_hline(yintercept =  0,   color = "black", linewidth = 0.5) +
    geom_hline(yintercept =  0.5, color = "black", linetype = "dashed", linewidth = 0.5) +
    geom_hline(yintercept = -0.5, color = "black", linetype = "dashed", linewidth = 0.5) +
    geom_point(size = 2, position = position_jitter(width = 0.1, seed = 42)) +
    scale_color_manual(values = colors_cell) +
    scale_y_continuous(limits = y_lims) +
    scale_x_discrete(drop = TRUE) +
    labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)"),
         title = paste0(gene_name, " (", protein_id, ")")) +
    theme_classic() +
    theme(
      plot.title      = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.title.y    = element_text(size = 9),
      axis.text.x     = element_text(color = "black", size = 9, angle = 90, hjust = 1),
      axis.text.y     = element_text(color = "black", size = 9),
      legend.position = "none",
      panel.spacing   = unit(0.3, "lines")
    )

  n_cells <- length(unique(df$cell))
  # base width = 1.1 in per cell + 0.6 in padding for axis
  w <- 0.6 + 1.1 * n_cells

  pdf_path <- paste0(out_dir, gene_name, '_', protein_id, '.pdf')
  png_path <- paste0(out_dir, gene_name, '_', protein_id, '.png')

  ggsave(pdf_path, plot = p, height = 2.5, width = w, units = 'in')
  ggsave(png_path, plot = p, height = 2.5, width = w, units = 'in', dpi = 600)

  cat(sprintf("Saved %s (%d cells, width=%.2f in)\n", gene_name, n_cells, w))
}

# ---- Generate all candidate plots -----------------------------------------
cat("\n=== Generating Figure 4C candidate slide plots ===\n")
for (i in seq_len(nrow(candidates))) {
  plot_candidate(candidates$Gene[i], candidates$Protein_ID[i])
}
cat("\nAll plots saved to:", out_dir, "\n")
