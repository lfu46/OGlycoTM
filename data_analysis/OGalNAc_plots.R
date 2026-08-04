# OGalNAc_plots.R
# 1. Site-level ranking plot (Intensity Ctrl/Tuni vs rank)
# 2. GSEA dotplot for exclusive O-GalNAc proteins

library(tidyverse)

source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'
figure_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/'

colors_cell <- c("HEK293T" = "#4DBBD5", "HepG2" = "#F39B7F", "Jurkat" = "#00A087")

dir.create(paste0(figure_file_path, "Figure6_OGalNAc"), showWarnings = FALSE)

# =============================================================================
# Part 1: Site-level ranking plot
# =============================================================================

# Aggregate PSMs to site level, compute ratio
all_sites <- data.frame()

for (cell in c("HEK293T", "HepG2", "Jurkat")) {
  site <- read_csv(
    paste0(source_file_path, "site/OGalNAc_site_", cell, ".csv"),
    show_col_types = FALSE
  )

  # Aggregate to unique site (sum intensities)
  site_quant <- site |>
    group_by(site_index) |>
    summarize(
      Gene = first(Gene),
      Protein.ID = first(Protein.ID),
      mean_Ctrl = mean(c(
        sum(Intensity.Ctrl_4), sum(Intensity.Ctrl_5), sum(Intensity.Ctrl_6)
      ) / 3),
      mean_Tuni = mean(c(
        sum(Intensity.Tuni_1), sum(Intensity.Tuni_2), sum(Intensity.Tuni_3)
      ) / 3),
      .groups = "drop"
    ) |>
    filter(mean_Ctrl > 0 & mean_Tuni > 0) |>
    mutate(
      ratio = mean_Ctrl / mean_Tuni,
      cell = cell
    ) |>
    arrange(ratio) |>
    mutate(rank = row_number())

  all_sites <- bind_rows(all_sites, site_quant)
  cat(cell, ":", nrow(site_quant), "quantified sites\n")
}

# Plot: each cell type as a separate panel
ranking_plot <- all_sites |>
  ggplot(aes(x = rank, y = ratio)) +
  geom_line(linewidth = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  facet_wrap(~ cell, scales = "free_x", nrow = 1) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    x = "Glycopeptide rank",
    y = "Intensity(Ctrl/Tuni)"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_blank()
  )

ggsave(
  paste0(figure_file_path, "Figure6_OGalNAc/OGalNAc_site_ranking.pdf"),
  ranking_plot, width = 7, height = 2.5, units = "in"
)
cat("Site ranking plot saved\n")

# =============================================================================
# Part 2: GSEA dotplot for exclusive O-GalNAc
# =============================================================================

# Combine significant GSEA results across cell types and ontologies
gsea_all <- data.frame()

for (cell in c("HEK293T", "HepG2", "Jurkat")) {
  for (ont in c("BP", "CC", "MF")) {
    path <- paste0(source_file_path, "enrichment/OGalNAc_exclusive_GSEA_", ont, "_", cell, ".csv")
    if (file.exists(path)) {
      df <- read_csv(path, show_col_types = FALSE)
      # Skip files with no significant results
      if (!"setSize" %in% names(df) || !is.numeric(df$setSize)) {
        df$setSize <- as.numeric(df$setSize)
      }
      df <- df |>
        filter(p.adjust < 0.05) |>
        mutate(cell = cell, ont = ont)
      if (nrow(df) > 0) {
        gsea_all <- bind_rows(gsea_all, df)
      }
    }
  }
}

cat("\nTotal significant GSEA terms:", nrow(gsea_all), "\n")

if (nrow(gsea_all) > 0) {
  # Jurkat only, both up and down
  selected_terms <- gsea_all |>
    filter(cell == "Jurkat") |>
    group_by(ont) |>
    slice_min(p.adjust, n = 4) |>
    ungroup() |>
    filter(!Description %in% c(
      "intracellular anatomical structure",
      "intracellular organelle",
      "organic cyclic compound binding"
    ))

  cat("Selected terms for plot:", nrow(selected_terms), "\n")

  # Create dotplot
  gsea_dotplot <- selected_terms |>
    mutate(
      Direction = ifelse(NES > 0, "Upregulated", "Downregulated"),
      Description = str_wrap(Description, width = 35),
      Description = reorder(Description, NES)
    ) |>
    ggplot(aes(x = NES, y = Description, size = setSize, color = p.adjust)) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    scale_color_gradient(low = "#004D40", high = "#80CBC4", name = "Adj. P") +
    scale_size_continuous(range = c(1.5, 6), name = "Gene set\nsize") +
    labs(
      x = "Normalized Enrichment Score (NES)",
      y = ""
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 9, color = "black"),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.text.x = element_text(size = 9, color = "black"),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.3, "cm"),
      panel.grid.minor = element_blank()
    )

  ggsave(
    paste0(figure_file_path, "Figure6_OGalNAc/OGalNAc_exclusive_GSEA_dotplot.pdf"),
    gsea_dotplot, width = 3.5, height = 2, units = "in"
  )
  cat("GSEA dotplot saved\n")
}

cat("\nPlots complete!\n")
