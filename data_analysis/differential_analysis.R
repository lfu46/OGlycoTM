# Differential Analysis for O-GlcNAc Protein and Site Level
# Using limma package for moderated t-test

library(tidyverse)
library(limma)

# Source data paths
source('data_source.R')

# =============================================================================
# Experimental Design
# =============================================================================

# Design matrix: Tuni (case1) vs Ctrl (case2), 3 replicates each
Experiment_Model <- model.matrix(
  ~ 0 + factor(rep(c("Tuni", "Ctrl"), each = 3), levels = c("Tuni", "Ctrl"))
)
colnames(Experiment_Model) <- c("Tuni", "Ctrl")

# Contrast: Tuni - Ctrl (positive logFC = higher in Tuni)
Contrast_Matrix <- makeContrasts(Tuni_vs_Ctrl = Tuni - Ctrl, levels = Experiment_Model)

# Define intensity columns (TMM normalized)
intensity_cols <- c(
  "Intensity.Tuni_1_sl_tmm", "Intensity.Tuni_2_sl_tmm", "Intensity.Tuni_3_sl_tmm",
  "Intensity.Ctrl_4_sl_tmm", "Intensity.Ctrl_5_sl_tmm", "Intensity.Ctrl_6_sl_tmm"
)

# =============================================================================
# Reusable Function for Limma Differential Analysis
# =============================================================================

run_limma_DE <- function(data, id_col, metadata_cols, intensity_cols,
                         design_matrix, contrast_matrix) {

  # Log2 transform intensity data
  log2_data <- data %>%
    mutate(across(all_of(intensity_cols), ~ log2(.x)))

  # Create data matrix for limma
  data_matrix <- log2_data %>%
    select(all_of(intensity_cols)) %>%
    as.matrix()

  # Set row names as identifiers
  rownames(data_matrix) <- data[[id_col]]

  # Fit linear model
  fit <- lmFit(data_matrix, design_matrix)

  # Apply contrasts
  fit_contrast <- contrasts.fit(fit, contrast_matrix)

  # Empirical Bayes moderation
  fit_contrast <- eBayes(fit_contrast)

  # Extract results
  top_table <- topTable(fit_contrast, number = Inf, adjust.method = "BH")

  # Convert to tibble and add identifier column
  result <- as_tibble(top_table)
  result[[id_col]] <- rownames(top_table)

  # Join with metadata
  metadata <- data %>%
    select(all_of(c(id_col, metadata_cols)))

  result <- result %>%
    left_join(metadata, by = id_col) %>%
    select(all_of(c(id_col, metadata_cols)), logFC, AveExpr, t, P.Value, adj.P.Val, B)

  return(result)
}

# =============================================================================
# Load Normalized Data
# =============================================================================

# Protein level
OGlcNAc_protein_norm_HEK293T <- read_csv(
  paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_HEK293T.csv')
)
OGlcNAc_protein_norm_HepG2 <- read_csv(
  paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_HepG2.csv')
)
OGlcNAc_protein_norm_Jurkat <- read_csv(
  paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_Jurkat.csv')
)

# Site level
OGlcNAc_site_norm_HEK293T <- read_csv(
  paste0(source_file_path, 'normalization/OGlcNAc_site_norm_HEK293T.csv')
)
OGlcNAc_site_norm_HepG2 <- read_csv(
  paste0(source_file_path, 'normalization/OGlcNAc_site_norm_HepG2.csv')
)
OGlcNAc_site_norm_Jurkat <- read_csv(
  paste0(source_file_path, 'normalization/OGlcNAc_site_norm_Jurkat.csv')
)

# =============================================================================
# Protein Level Differential Analysis
# =============================================================================

protein_metadata_cols <- c("Entry.Name", "Gene", "Protein.Description")

OGlcNAc_protein_DE_HEK293T <- run_limma_DE(
  data = OGlcNAc_protein_norm_HEK293T,
  id_col = "Protein.ID",
  metadata_cols = protein_metadata_cols,
  intensity_cols = intensity_cols,
  design_matrix = Experiment_Model,
  contrast_matrix = Contrast_Matrix
)

OGlcNAc_protein_DE_HepG2 <- run_limma_DE(
  data = OGlcNAc_protein_norm_HepG2,
  id_col = "Protein.ID",
  metadata_cols = protein_metadata_cols,
  intensity_cols = intensity_cols,
  design_matrix = Experiment_Model,
  contrast_matrix = Contrast_Matrix
)

OGlcNAc_protein_DE_Jurkat <- run_limma_DE(
  data = OGlcNAc_protein_norm_Jurkat,
  id_col = "Protein.ID",
  metadata_cols = protein_metadata_cols,
  intensity_cols = intensity_cols,
  design_matrix = Experiment_Model,
  contrast_matrix = Contrast_Matrix
)

# =============================================================================
# Site Level Differential Analysis
# =============================================================================

site_metadata_cols <- c("Protein.ID", "Entry.Name", "Gene", "Protein.Description")

OGlcNAc_site_DE_HEK293T <- run_limma_DE(
  data = OGlcNAc_site_norm_HEK293T,
  id_col = "site_index",
  metadata_cols = site_metadata_cols,
  intensity_cols = intensity_cols,
  design_matrix = Experiment_Model,
  contrast_matrix = Contrast_Matrix
)

OGlcNAc_site_DE_HepG2 <- run_limma_DE(
  data = OGlcNAc_site_norm_HepG2,
  id_col = "site_index",
  metadata_cols = site_metadata_cols,
  intensity_cols = intensity_cols,
  design_matrix = Experiment_Model,
  contrast_matrix = Contrast_Matrix
)

OGlcNAc_site_DE_Jurkat <- run_limma_DE(
  data = OGlcNAc_site_norm_Jurkat,
  id_col = "site_index",
  metadata_cols = site_metadata_cols,
  intensity_cols = intensity_cols,
  design_matrix = Experiment_Model,
  contrast_matrix = Contrast_Matrix
)

# =============================================================================
# Save Results
# =============================================================================

# Create output directory if it doesn't exist
output_dir <- paste0(source_file_path, 'differential_analysis/')
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save protein level results
write_csv(OGlcNAc_protein_DE_HEK293T, paste0(output_dir, 'OGlcNAc_protein_DE_HEK293T.csv'))
write_csv(OGlcNAc_protein_DE_HepG2, paste0(output_dir, 'OGlcNAc_protein_DE_HepG2.csv'))
write_csv(OGlcNAc_protein_DE_Jurkat, paste0(output_dir, 'OGlcNAc_protein_DE_Jurkat.csv'))

# Save site level results
write_csv(OGlcNAc_site_DE_HEK293T, paste0(output_dir, 'OGlcNAc_site_DE_HEK293T.csv'))
write_csv(OGlcNAc_site_DE_HepG2, paste0(output_dir, 'OGlcNAc_site_DE_HepG2.csv'))
write_csv(OGlcNAc_site_DE_Jurkat, paste0(output_dir, 'OGlcNAc_site_DE_Jurkat.csv'))

cat("Differential analysis complete!\n")
cat("Results saved to:", output_dir, "\n")

# =============================================================================
# Summary Statistics
# =============================================================================

summarize_DE <- function(de_result, name, fc_threshold = 0.5, pval_threshold = 0.05) {
  sig_up <- de_result %>% filter(logFC > fc_threshold & adj.P.Val < pval_threshold) %>% nrow()
  sig_down <- de_result %>% filter(logFC < -fc_threshold & adj.P.Val < pval_threshold) %>% nrow()
  total <- nrow(de_result)

  cat(sprintf("%s: %d total, %d up, %d down (|logFC| > %g, adj.P.Val < %g)\n",
              name, total, sig_up, sig_down, fc_threshold, pval_threshold))
}

cat("\n=== Summary (Tuni vs Ctrl) ===\n")
cat("Protein Level:\n")
summarize_DE(OGlcNAc_protein_DE_HEK293T, "  HEK293T")
summarize_DE(OGlcNAc_protein_DE_HepG2, "  HepG2")
summarize_DE(OGlcNAc_protein_DE_Jurkat, "  Jurkat")

cat("\nSite Level:\n")
summarize_DE(OGlcNAc_site_DE_HEK293T, "  HEK293T")
summarize_DE(OGlcNAc_site_DE_HepG2, "  HepG2")
summarize_DE(OGlcNAc_site_DE_Jurkat, "  Jurkat")

# =============================================================================
# Commonly Up/Downregulated Proteins Across Cell Types
# =============================================================================

# Function to classify proteins as "up", "down", or "ns" (not significant)
classify_DE <- function(de_result, fc_threshold = 0.5, pval_threshold = 0.05) {
  de_result %>%
    mutate(
      DE_status = case_when(
        logFC > fc_threshold & adj.P.Val < pval_threshold ~ "up",
        logFC < -fc_threshold & adj.P.Val < pval_threshold ~ "down",
        TRUE ~ "ns"
      )
    ) %>%
    select(Protein.ID, DE_status)
}

# Classify proteins in each cell type
DE_class_HEK293T <- classify_DE(OGlcNAc_protein_DE_HEK293T) %>%
  rename(DE_HEK293T = DE_status)
DE_class_HepG2 <- classify_DE(OGlcNAc_protein_DE_HepG2) %>%
  rename(DE_HepG2 = DE_status)
DE_class_Jurkat <- classify_DE(OGlcNAc_protein_DE_Jurkat) %>%
  rename(DE_Jurkat = DE_status)

# Combine classifications
DE_combined <- DE_class_HEK293T %>%
  full_join(DE_class_HepG2, by = "Protein.ID") %>%
  full_join(DE_class_Jurkat, by = "Protein.ID")

# Count up/down occurrences for each protein
DE_combined <- DE_combined %>%
  mutate(
    n_up = (DE_HEK293T == "up") + (DE_HepG2 == "up") + (DE_Jurkat == "up"),
    n_down = (DE_HEK293T == "down") + (DE_HepG2 == "down") + (DE_Jurkat == "down")
  )

# Commonly upregulated: up in ≥2 cell types AND down in 0 cell types
commonly_up <- DE_combined %>%
  filter(n_up >= 2 & n_down == 0)

# Commonly downregulated: down in ≥2 cell types AND up in 0 cell types
commonly_down <- DE_combined %>%
  filter(n_down >= 2 & n_up == 0)

cat("\n=== Commonly Regulated Proteins (|logFC| > 0.5, adj.P.Val < 0.05) ===\n")
cat("Definition: Up/down in ≥2 cell types AND not down/up in the other\n\n")
cat("Commonly upregulated proteins:", nrow(commonly_up), "\n")
cat("  - Up in all 3 cell types:", sum(commonly_up$n_up == 3), "\n")
cat("  - Up in 2 cell types:", sum(commonly_up$n_up == 2), "\n")
cat("\nCommonly downregulated proteins:", nrow(commonly_down), "\n")
cat("  - Down in all 3 cell types:", sum(commonly_down$n_down == 3), "\n")
cat("  - Down in 2 cell types:", sum(commonly_down$n_down == 2), "\n")

# Save commonly regulated protein lists
write_csv(commonly_up, paste0(output_dir, 'OGlcNAc_protein_commonly_up.csv'))
write_csv(commonly_down, paste0(output_dir, 'OGlcNAc_protein_commonly_down.csv'))

cat("\nCommonly regulated protein lists saved to:", output_dir, "\n")


# =============================================================================
# Whole Proteome Differential Analysis
# =============================================================================

# Load normalized WP data
WP_protein_norm_HEK293T <- read_csv(
  paste0(source_file_path, 'normalization/WP_protein_norm_HEK293T.csv')
)
WP_protein_norm_HepG2 <- read_csv(
  paste0(source_file_path, 'normalization/WP_protein_norm_HepG2.csv')
)
WP_protein_norm_Jurkat <- read_csv(
  paste0(source_file_path, 'normalization/WP_protein_norm_Jurkat.csv')
)

# Define metadata columns for WP
WP_metadata_cols <- c("Gene.Symbol", "Protein.ID", "Protein.MWT.kDa.", "Annotation")

# HEK293T
WP_protein_DE_HEK293T <- run_limma_DE(
  data = WP_protein_norm_HEK293T,
  id_col = "UniProt_Accession",
  metadata_cols = WP_metadata_cols,
  intensity_cols = intensity_cols,
  design_matrix = Experiment_Model,
  contrast_matrix = Contrast_Matrix
)

# HepG2
WP_protein_DE_HepG2 <- run_limma_DE(
  data = WP_protein_norm_HepG2,
  id_col = "UniProt_Accession",
  metadata_cols = WP_metadata_cols,
  intensity_cols = intensity_cols,
  design_matrix = Experiment_Model,
  contrast_matrix = Contrast_Matrix
)

# Jurkat
WP_protein_DE_Jurkat <- run_limma_DE(
  data = WP_protein_norm_Jurkat,
  id_col = "UniProt_Accession",
  metadata_cols = WP_metadata_cols,
  intensity_cols = intensity_cols,
  design_matrix = Experiment_Model,
  contrast_matrix = Contrast_Matrix
)

# Save WP differential analysis results
write_csv(WP_protein_DE_HEK293T, paste0(source_file_path, 'differential_analysis/WP_protein_DE_HEK293T.csv'))
write_csv(WP_protein_DE_HepG2, paste0(source_file_path, 'differential_analysis/WP_protein_DE_HepG2.csv'))
write_csv(WP_protein_DE_Jurkat, paste0(source_file_path, 'differential_analysis/WP_protein_DE_Jurkat.csv'))

cat("\nWhole Proteome differential analysis complete!\n")

# =============================================================================
# Whole Proteome Summary Statistics
# =============================================================================

cat("\n=== Whole Proteome Summary (Tuni vs Ctrl) ===\n")
summarize_DE(WP_protein_DE_HEK293T, "  HEK293T")
summarize_DE(WP_protein_DE_HepG2, "  HepG2")
summarize_DE(WP_protein_DE_Jurkat, "  Jurkat")

# =============================================================================
# Whole Proteome Commonly Regulated Proteins
# =============================================================================

# Function to classify WP proteins (uses UniProt_Accession as ID)
classify_WP_DE <- function(de_result, fc_threshold = 0.5, pval_threshold = 0.05) {
  de_result %>%
    mutate(
      DE_status = case_when(
        logFC > fc_threshold & adj.P.Val < pval_threshold ~ "up",
        logFC < -fc_threshold & adj.P.Val < pval_threshold ~ "down",
        TRUE ~ "ns"
      )
    ) %>%
    select(UniProt_Accession, DE_status)
}

# Classify WP proteins in each cell type
WP_DE_class_HEK293T <- classify_WP_DE(WP_protein_DE_HEK293T) %>%
  rename(DE_HEK293T = DE_status)
WP_DE_class_HepG2 <- classify_WP_DE(WP_protein_DE_HepG2) %>%
  rename(DE_HepG2 = DE_status)
WP_DE_class_Jurkat <- classify_WP_DE(WP_protein_DE_Jurkat) %>%
  rename(DE_Jurkat = DE_status)

# Combine classifications
WP_DE_combined <- WP_DE_class_HEK293T %>%
  full_join(WP_DE_class_HepG2, by = "UniProt_Accession") %>%
  full_join(WP_DE_class_Jurkat, by = "UniProt_Accession")

# Count up/down occurrences for each protein
WP_DE_combined <- WP_DE_combined %>%
  mutate(
    n_up = (DE_HEK293T == "up") + (DE_HepG2 == "up") + (DE_Jurkat == "up"),
    n_down = (DE_HEK293T == "down") + (DE_HepG2 == "down") + (DE_Jurkat == "down")
  )

# Commonly upregulated: up in ≥2 cell types AND down in 0 cell types
WP_commonly_up <- WP_DE_combined %>%
  filter(n_up >= 2 & n_down == 0)

# Commonly downregulated: down in ≥2 cell types AND up in 0 cell types
WP_commonly_down <- WP_DE_combined %>%
  filter(n_down >= 2 & n_up == 0)

cat("\n=== Whole Proteome Commonly Regulated Proteins (|logFC| > 0.5, adj.P.Val < 0.05) ===\n")
cat("Definition: Up/down in ≥2 cell types AND not down/up in the other\n\n")
cat("Commonly upregulated proteins:", nrow(WP_commonly_up), "\n")
cat("  - Up in all 3 cell types:", sum(WP_commonly_up$n_up == 3), "\n")
cat("  - Up in 2 cell types:", sum(WP_commonly_up$n_up == 2), "\n")
cat("\nCommonly downregulated proteins:", nrow(WP_commonly_down), "\n")
cat("  - Down in all 3 cell types:", sum(WP_commonly_down$n_down == 3), "\n")
cat("  - Down in 2 cell types:", sum(WP_commonly_down$n_down == 2), "\n")

# Save WP commonly regulated protein lists
write_csv(WP_commonly_up, paste0(source_file_path, 'differential_analysis/WP_protein_commonly_up.csv'))
write_csv(WP_commonly_down, paste0(source_file_path, 'differential_analysis/WP_protein_commonly_down.csv'))

cat("\nWhole Proteome commonly regulated protein lists saved to:", source_file_path, "differential_analysis/\n")
