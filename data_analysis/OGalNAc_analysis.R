# OGalNAc_analysis.R
# Complete O-GalNAc analysis pipeline: quantification → normalization → DE → GO
# Following the same pipeline as O-GlcNAc

library(tidyverse)
library(edgeR)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)

# =============================================================================
# Source data
# =============================================================================
source("data_source.R")

source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'
figure_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/'

intensity_cols <- c("Intensity.Tuni_1", "Intensity.Tuni_2", "Intensity.Tuni_3",
                    "Intensity.Ctrl_4", "Intensity.Ctrl_5", "Intensity.Ctrl_6")

OGalNAc_compositions <- c(
  'HexNAt(1)GAO_Methoxylamine(1) % 326.1339',
  'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968'
)

# =============================================================================
# Step 1: Filter bonafide for O-GalNAc
# =============================================================================
cat(strrep("=", 60), "\nStep 1: O-GalNAc Filtering\n", strrep("=", 60), "\n")

OGalNAc_HEK293T <- OGlyco_HEK293T_bonafide |>
  filter(Total.Glycan.Composition %in% OGalNAc_compositions)
OGalNAc_HepG2 <- OGlyco_HepG2_bonafide |>
  filter(Total.Glycan.Composition %in% OGalNAc_compositions)
OGalNAc_Jurkat <- OGlyco_Jurkat_bonafide |>
  filter(Total.Glycan.Composition %in% OGalNAc_compositions)

write_csv(OGalNAc_HEK293T, paste0(source_file_path, "filtered/OGalNAc_HEK293T.csv"))
write_csv(OGalNAc_HepG2, paste0(source_file_path, "filtered/OGalNAc_HepG2.csv"))
write_csv(OGalNAc_Jurkat, paste0(source_file_path, "filtered/OGalNAc_Jurkat.csv"))

for (cell in c("HEK293T", "HepG2", "Jurkat")) {
  d <- get(paste0("OGalNAc_", cell))
  cat(cell, ": ", nrow(d), " PSMs, ", n_distinct(d$Protein.ID), " proteins\n")
}

# =============================================================================
# Step 2: Protein-level quantification (sum PSM intensities per protein)
# =============================================================================
cat("\n", strrep("=", 60), "\nStep 2: Protein-Level Quantification\n", strrep("=", 60), "\n")

quantify_protein <- function(data) {
  # Get metadata (first row per protein)
  meta <- data |>
    group_by(Protein.ID) |>
    slice_head(n = 1) |>
    ungroup() |>
    dplyr::select(Protein.ID, Entry.Name, Gene, Protein.Description)

  # Sum intensities
  quant <- data |>
    group_by(Protein.ID) |>
    summarize(
      across(starts_with("Intensity."), ~ sum(.x, na.rm = TRUE)),
      .groups = "drop"
    )

  left_join(quant, meta, by = "Protein.ID")
}

OGalNAc_protein_quant_HEK293T <- quantify_protein(OGalNAc_HEK293T)
OGalNAc_protein_quant_HepG2 <- quantify_protein(OGalNAc_HepG2)
OGalNAc_protein_quant_Jurkat <- quantify_protein(OGalNAc_Jurkat)

write_csv(OGalNAc_protein_quant_HEK293T, paste0(source_file_path, "quantification/OGalNAc_protein_quant_HEK293T.csv"))
write_csv(OGalNAc_protein_quant_HepG2, paste0(source_file_path, "quantification/OGalNAc_protein_quant_HepG2.csv"))
write_csv(OGalNAc_protein_quant_Jurkat, paste0(source_file_path, "quantification/OGalNAc_protein_quant_Jurkat.csv"))

for (cell in c("HEK293T", "HepG2", "Jurkat")) {
  cat(cell, ": ", nrow(get(paste0("OGalNAc_protein_quant_", cell))), " proteins quantified\n")
}

# =============================================================================
# Step 3: Normalization (SL + TMM)
# =============================================================================
cat("\n", strrep("=", 60), "\nStep 3: Normalization\n", strrep("=", 60), "\n")

filter_zero_na <- function(data, intensity_cols) {
  data |> filter(if_all(all_of(intensity_cols), ~ . > 0 & !is.na(.)))
}

normalize_intensity <- function(data, intensity_cols) {
  intensity_matrix <- data |> dplyr::select(all_of(intensity_cols)) |> as.matrix()

  # SL normalization
  col_sums <- colSums(intensity_matrix)
  target_mean <- mean(col_sums)
  sl_factors <- target_mean / col_sums
  intensity_sl <- sweep(intensity_matrix, 2, sl_factors, FUN = "*")
  intensity_sl_tb <- as_tibble(intensity_sl)
  colnames(intensity_sl_tb) <- paste0(intensity_cols, "_sl")

  # TMM normalization
  tmm_factors <- calcNormFactors(intensity_sl)
  intensity_sl_tmm <- sweep(intensity_sl, 2, tmm_factors, FUN = "/")
  intensity_sl_tmm_tb <- as_tibble(intensity_sl_tmm)
  colnames(intensity_sl_tmm_tb) <- paste0(intensity_cols, "_sl_tmm")

  bind_cols(data, intensity_sl_tb, intensity_sl_tmm_tb)
}

for (cell in c("HEK293T", "HepG2", "Jurkat")) {
  quant <- get(paste0("OGalNAc_protein_quant_", cell))
  filtered <- filter_zero_na(quant, intensity_cols)
  normed <- normalize_intensity(filtered, intensity_cols)
  assign(paste0("OGalNAc_protein_norm_", cell), normed)
  write_csv(normed, paste0(source_file_path, "normalization/OGalNAc_protein_norm_", cell, ".csv"))
  cat(cell, ": ", nrow(quant), " -> ", nrow(filtered), " after filtering -> normalized\n")
}

# =============================================================================
# Step 4: Differential Analysis (limma)
# =============================================================================
cat("\n", strrep("=", 60), "\nStep 4: Differential Analysis\n", strrep("=", 60), "\n")

Experiment_Model <- model.matrix(
  ~ 0 + factor(rep(c("Tuni", "Ctrl"), each = 3), levels = c("Tuni", "Ctrl"))
)
colnames(Experiment_Model) <- c("Tuni", "Ctrl")
Contrast_Matrix <- makeContrasts(Tuni_vs_Ctrl = Tuni - Ctrl, levels = Experiment_Model)

norm_intensity_cols <- paste0(intensity_cols, "_sl_tmm")

run_limma_DE <- function(data, id_col, metadata_cols, intensity_cols, design_matrix, contrast_matrix) {
  log2_data <- data %>% mutate(across(all_of(intensity_cols), ~ log2(.x)))
  data_matrix <- log2_data %>% dplyr::select(all_of(intensity_cols)) %>% as.matrix()
  rownames(data_matrix) <- data[[id_col]]
  fit <- lmFit(data_matrix, design_matrix)
  fit_contrast <- contrasts.fit(fit, contrast_matrix)
  fit_contrast <- eBayes(fit_contrast)
  top_table <- topTable(fit_contrast, number = Inf, adjust.method = "BH")
  result <- as_tibble(top_table)
  result[[id_col]] <- rownames(top_table)
  metadata <- data %>% dplyr::select(all_of(c(id_col, metadata_cols)))
  result %>%
    left_join(metadata, by = id_col) %>%
    dplyr::select(all_of(c(id_col, metadata_cols)), logFC, AveExpr, t, P.Value, adj.P.Val, B)
}

for (cell in c("HEK293T", "HepG2", "Jurkat")) {
  normed <- get(paste0("OGalNAc_protein_norm_", cell))
  de <- run_limma_DE(normed, "Protein.ID", c("Gene", "Protein.Description"),
                     norm_intensity_cols, Experiment_Model, Contrast_Matrix)
  assign(paste0("OGalNAc_protein_DE_", cell), de)
  write_csv(de, paste0(source_file_path, "differential_analysis/OGalNAc_protein_DE_", cell, ".csv"))

  n_up <- sum(de$logFC > 0.5 & de$adj.P.Val < 0.05)
  n_down <- sum(de$logFC < -0.5 & de$adj.P.Val < 0.05)
  cat(cell, ": ", nrow(de), " proteins, ", n_up, " up, ", n_down, " down\n")
}

# =============================================================================
# Step 5: GO Enrichment of regulated O-GalNAc proteins
# =============================================================================
cat("\n", strrep("=", 60), "\nStep 5: GO Enrichment\n", strrep("=", 60), "\n")

for (cell in c("HEK293T", "HepG2", "Jurkat")) {
  de <- get(paste0("OGalNAc_protein_DE_", cell))

  sig_up <- de |> filter(logFC > 0.5, adj.P.Val < 0.05) |> pull(Protein.ID)
  sig_down <- de |> filter(logFC < -0.5, adj.P.Val < 0.05) |> pull(Protein.ID)
  universe <- de$Protein.ID

  cat("\n", cell, ":\n")
  cat("  Upregulated: ", length(sig_up), "\n")
  cat("  Downregulated: ", length(sig_down), "\n")

  # Upregulated GO
  if (length(sig_up) >= 5) {
    go_up <- enrichGO(
      gene = sig_up, OrgDb = org.Hs.eg.db, universe = universe,
      keyType = 'UNIPROT', ont = 'ALL', pvalueCutoff = 0.05, qvalueCutoff = 0.1
    )
    write_csv(go_up@result, paste0(source_file_path, "enrichment/OGalNAc_up_GO_", cell, ".csv"))
    cat("  Upregulated GO terms (p<0.05): ", sum(go_up@result$pvalue < 0.05), "\n")
  }

  # Downregulated GO
  if (length(sig_down) >= 5) {
    go_down <- enrichGO(
      gene = sig_down, OrgDb = org.Hs.eg.db, universe = universe,
      keyType = 'UNIPROT', ont = 'ALL', pvalueCutoff = 0.05, qvalueCutoff = 0.1
    )
    write_csv(go_down@result, paste0(source_file_path, "enrichment/OGalNAc_down_GO_", cell, ".csv"))
    cat("  Downregulated GO terms (p<0.05): ", sum(go_down@result$pvalue < 0.05), "\n")
  }
}

cat("\n", strrep("=", 60), "\nO-GalNAc analysis complete!\n", strrep("=", 60), "\n")
