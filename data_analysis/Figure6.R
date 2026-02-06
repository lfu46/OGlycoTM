# Figure 6: O-GlcNAc Site Structural Feature Analysis
#
# This analysis investigates structural and sequence features associated with
# tunicamycin-induced O-GlcNAc fold changes across three cell types.
#
# Key Finding: O-GlcNAc sites in structured regions (pLDDT >= 50) show larger
# treatment responses compared to sites in IDR regions (pLDDT < 50).
#
# Figure Panels:
#   - Figure 6A: logFC distribution by cell type
#   - Figure 6B: Secondary structure distribution (donut plot)
#   - Figure 6C: Standardized regression coefficients
#   - Figure 6D: IDR vs Structured effect across cell types
#   - Figure 6E: HYOU1 T871 - Structured region example (logFC = +2.47, pLDDT = 94.1)
#   - Figure 6F: HOXA13 S199 - IDR region example (logFC = -0.02, pLDDT = 41.6)
#
# Data files generated:
#   - site_features/OGlcNAc_site_features.csv
#   - site_features/regression_coefficients_standardized.csv
#
# Data files required:
#   - site_features/alphafold_structural_features.csv (pLDDT profiles)

# ==============================================================================
# PART 1: DATA PREPARATION AND FEATURE EXTRACTION
# ==============================================================================

# Load data sources
source('data_source.R')
source('data_source_DE.R')

# Additional packages
library(Peptides)
library(ggpubr)

# ------------------------------------------------------------------------------
# Section 1.1: Load Site-Level Differential Analysis Data
# ------------------------------------------------------------------------------

OGlcNAc_site_DE_HEK293T <- read_csv(
  paste0(source_file_path, 'differential_analysis/OGlcNAc_site_DE_HEK293T.csv'),
  show_col_types = FALSE
)
OGlcNAc_site_DE_HepG2 <- read_csv(
  paste0(source_file_path, 'differential_analysis/OGlcNAc_site_DE_HepG2.csv'),
  show_col_types = FALSE
)
OGlcNAc_site_DE_Jurkat <- read_csv(
  paste0(source_file_path, 'differential_analysis/OGlcNAc_site_DE_Jurkat.csv'),
  show_col_types = FALSE
)

cat("=== Site-Level Differential Expression Data Summary ===\n")
cat("HEK293T:", nrow(OGlcNAc_site_DE_HEK293T), "sites\n")
cat("HepG2:", nrow(OGlcNAc_site_DE_HepG2), "sites\n")
cat("Jurkat:", nrow(OGlcNAc_site_DE_Jurkat), "sites\n\n")

# ------------------------------------------------------------------------------
# Section 1.2: Extract Basic Features
# ------------------------------------------------------------------------------

extract_basic_features <- function(site_data, de_data) {
  site_info <- site_data |>
    group_by(site_index) |>
    slice_head(n = 1) |>
    ungroup() |>
    dplyr::select(
      site_index, Protein.ID, Gene, Peptide,
      Protein.Start, Protein.End, modified_residue, site_number
    )

  features <- de_data |>
    left_join(site_info, by = c("site_index", "Protein.ID", "Gene")) |>
    mutate(is_serine = as.integer(modified_residue == "S"))

  return(features)
}

features_HEK293T <- extract_basic_features(OGlcNAc_site_HEK293T, OGlcNAc_site_DE_HEK293T)
features_HepG2 <- extract_basic_features(OGlcNAc_site_HepG2, OGlcNAc_site_DE_HepG2)
features_Jurkat <- extract_basic_features(OGlcNAc_site_Jurkat, OGlcNAc_site_DE_Jurkat)

# ------------------------------------------------------------------------------
# Section 1.3: Calculate sites_per_protein
# ------------------------------------------------------------------------------

add_sites_per_protein <- function(features_df) {
  sites_count <- features_df |>
    group_by(Protein.ID) |>
    summarize(sites_per_protein = n(), .groups = "drop")

  features_df |>
    left_join(sites_count, by = "Protein.ID")
}

features_HEK293T <- add_sites_per_protein(features_HEK293T)
features_HepG2 <- add_sites_per_protein(features_HepG2)
features_Jurkat <- add_sites_per_protein(features_Jurkat)

# ------------------------------------------------------------------------------
# Section 1.4: Extract 7-mer Sequence and Calculate Sequence Properties
# ------------------------------------------------------------------------------

library(Biostrings)

fasta_file <- paste0(source_file_path, 'reference/uniprotkb_reviewed_true_AND_model_organ_2026_01_09.fasta')
protein_sequences <- readAAStringSet(fasta_file)

parse_uniprot_id <- function(header) {
  parts <- strsplit(header, "\\|")[[1]]
  if (length(parts) >= 2) return(parts[2])
  return(NA_character_)
}

protein_ids <- sapply(names(protein_sequences), parse_uniprot_id)
names(protein_sequences) <- protein_ids
protein_seq_lookup <- as.character(protein_sequences)
names(protein_seq_lookup) <- protein_ids

cat("Loaded", length(protein_seq_lookup), "protein sequences from FASTA\n")

extract_7mer_from_protein <- function(protein_id, site_number, seq_lookup) {
  if (is.na(protein_id) | is.na(site_number)) return(NA_character_)
  protein_seq <- seq_lookup[protein_id]
  if (is.na(protein_seq)) return(NA_character_)

  seq_len <- nchar(protein_seq)
  site_pos <- as.integer(site_number)
  if (site_pos < 1 | site_pos > seq_len) return(NA_character_)

  start_pos <- max(1, site_pos - 3)
  end_pos <- min(seq_len, site_pos + 3)
  seq_7mer <- substr(protein_seq, start_pos, end_pos)

  if (site_pos <= 3) {
    left_pad <- paste(rep("_", 3 - site_pos + 1), collapse = "")
    seq_7mer <- paste0(left_pad, seq_7mer)
  }
  if (site_pos > seq_len - 3) {
    right_pad <- paste(rep("_", 3 - (seq_len - site_pos)), collapse = "")
    seq_7mer <- paste0(seq_7mer, right_pad)
  }
  return(seq_7mer)
}

calculate_sequence_properties <- function(seq_7mer) {
  if (is.na(seq_7mer) | nchar(seq_7mer) < 3) {
    return(list(pI = NA_real_, hydrophobicity = NA_real_))
  }
  clean_seq <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", toupper(seq_7mer))
  if (nchar(clean_seq) < 3) {
    return(list(pI = NA_real_, hydrophobicity = NA_real_))
  }
  pI_val <- tryCatch(pI(clean_seq, pKscale = "EMBOSS"), error = function(e) NA_real_)
  hydro_val <- tryCatch(hydrophobicity(clean_seq, scale = "KyteDoolittle"), error = function(e) NA_real_)
  return(list(pI = pI_val, hydrophobicity = hydro_val))
}

add_sequence_features <- function(features_df, seq_lookup) {
  features_df <- features_df |>
    rowwise() |>
    mutate(seq_7mer = extract_7mer_from_protein(Protein.ID, site_number, seq_lookup)) |>
    ungroup()

  features_df <- features_df |>
    rowwise() |>
    mutate(seq_props = list(calculate_sequence_properties(seq_7mer))) |>
    ungroup() |>
    mutate(
      pI_7mer = sapply(seq_props, function(x) x$pI),
      hydrophobicity_7mer = sapply(seq_props, function(x) x$hydrophobicity)
    ) |>
    dplyr::select(-seq_props)

  return(features_df)
}

cat("Extracting 7-mer sequences and calculating properties...\n")
features_HEK293T <- add_sequence_features(features_HEK293T, protein_seq_lookup)
features_HepG2 <- add_sequence_features(features_HepG2, protein_seq_lookup)
features_Jurkat <- add_sequence_features(features_Jurkat, protein_seq_lookup)

# ------------------------------------------------------------------------------
# Section 1.5: Add Protein Abundance from Whole Proteome
# ------------------------------------------------------------------------------

add_protein_abundance <- function(features_df, wp_data) {
  abundance_data <- wp_data[, c("UniProt_Accession", "AveExpr")]
  names(abundance_data)[2] <- "log_protein_abundance"
  features_df |>
    left_join(abundance_data, by = c("Protein.ID" = "UniProt_Accession"))
}

features_HEK293T <- add_protein_abundance(features_HEK293T, WP_protein_DE_HEK293T)
features_HepG2 <- add_protein_abundance(features_HepG2, WP_protein_DE_HepG2)
features_Jurkat <- add_protein_abundance(features_Jurkat, WP_protein_DE_Jurkat)

# ------------------------------------------------------------------------------
# Section 1.6: Extract Structural Features from AlphaFold (StructureMap)
# ------------------------------------------------------------------------------

library(reticulate)

alphafold_cif_folder <- paste0(source_file_path, "reference/alphafold_cif/")
alphafold_pae_folder <- paste0(source_file_path, "reference/alphafold_pae/")
structuremap_env_name <- "structuremap_env"

env_exists <- tryCatch({
  conda_list() |> filter(name == structuremap_env_name) |> nrow() > 0
}, error = function(e) FALSE)

if (!env_exists) {
  cat("=== StructureMap Section SKIPPED ===\n")
  cat("Conda environment '", structuremap_env_name, "' not found.\n", sep = "")
  features_HEK293T <- features_HEK293T |>
    mutate(pPSE_24 = NA_real_, pPSE_12 = NA_real_, pPSE_24_smooth10 = NA_real_,
           is_IDR = NA_integer_, secondary_structure_simple = NA_character_, pLDDT = NA_real_)
  features_HepG2 <- features_HepG2 |>
    mutate(pPSE_24 = NA_real_, pPSE_12 = NA_real_, pPSE_24_smooth10 = NA_real_,
           is_IDR = NA_integer_, secondary_structure_simple = NA_character_, pLDDT = NA_real_)
  features_Jurkat <- features_Jurkat |>
    mutate(pPSE_24 = NA_real_, pPSE_12 = NA_real_, pPSE_24_smooth10 = NA_real_,
           is_IDR = NA_integer_, secondary_structure_simple = NA_character_, pLDDT = NA_real_)
} else {
  use_condaenv(structuremap_env_name, required = TRUE)
  source_python("structuremap_analysis.py")

  all_proteins <- unique(c(
    features_HEK293T$Protein.ID,
    features_HepG2$Protein.ID,
    features_Jurkat$Protein.ID
  ))
  all_proteins <- all_proteins[!is.na(all_proteins)]

  cat("=== StructureMap Structural Feature Extraction ===\n")
  cat("Total unique proteins:", length(all_proteins), "\n\n")

  structural_cache_file <- paste0(source_file_path, 'site_features/alphafold_structural_features.csv')

  if (file.exists(structural_cache_file)) {
    cat("Loading cached structural features from:", structural_cache_file, "\n")
    structural_features <- read_csv(structural_cache_file, show_col_types = FALSE)
  } else {
    cat("Step 1: SKIPPED - AlphaFold files already downloaded\n\n")
    cat("Step 2: Extracting structural features...\n")
    structural_features <- extract_structural_features(
      protein_ids = all_proteins,
      cif_folder = alphafold_cif_folder,
      pae_folder = alphafold_pae_folder
    )
    structural_features <- as_tibble(structural_features)
    write_csv(structural_features, structural_cache_file)
    cat("\nStructural features saved to:", structural_cache_file, "\n")
  }

  merge_structural_features <- function(features_df, structural_df) {
    features_df |>
      left_join(
        structural_df |>
          dplyr::select(protein_id, position, pLDDT, pPSE_24, pPSE_12,
                        pPSE_24_smooth10, is_IDR, secondary_structure_simple),
        by = c("Protein.ID" = "protein_id", "site_number" = "position")
      )
  }

  features_HEK293T <- merge_structural_features(features_HEK293T, structural_features)
  features_HepG2 <- merge_structural_features(features_HepG2, structural_features)
  features_Jurkat <- merge_structural_features(features_Jurkat, structural_features)

  cat("\n=== Structural Feature Coverage ===\n")
  cat("HEK293T:", sum(!is.na(features_HEK293T$pPSE_24)), "/", nrow(features_HEK293T), "sites\n")
  cat("HepG2:", sum(!is.na(features_HepG2$pPSE_24)), "/", nrow(features_HepG2), "sites\n")
  cat("Jurkat:", sum(!is.na(features_Jurkat$pPSE_24)), "/", nrow(features_Jurkat), "sites\n")
}

# ------------------------------------------------------------------------------
# Section 1.7: Combine All Features into Analysis Dataframe
# ------------------------------------------------------------------------------

features_combined <- bind_rows(
  features_HEK293T |> mutate(cell = "HEK293T"),
  features_HepG2 |> mutate(cell = "HepG2"),
  features_Jurkat |> mutate(cell = "Jurkat")
)

features_analysis <- features_combined |>
  dplyr::select(
    site_index, Protein.ID, Gene, cell, site_number, Peptide, seq_7mer,
    logFC, adj.P.Val,
    is_serine, sites_per_protein, pI_7mer, hydrophobicity_7mer, log_protein_abundance,
    pLDDT, pPSE_24, pPSE_12, pPSE_24_smooth10, is_IDR, secondary_structure_simple
  )

cat("\n=== Final Analysis Dataframe Summary ===\n")
cat("Total observations:", nrow(features_analysis), "\n")
cat("Complete cases:", sum(complete.cases(features_analysis)), "\n")

write_csv(
  features_analysis,
  paste0(source_file_path, 'site_features/OGlcNAc_site_features.csv')
)
cat("Features saved to:", paste0(source_file_path, 'site_features/OGlcNAc_site_features.csv'), "\n")


# ==============================================================================
# PART 2: REGRESSION ANALYSIS
# ==============================================================================

cat("\n=== Combined Regression with Cell Type Factor ===\n\n")

# Prepare combined regression data
reg_data_combined <- features_analysis |>
  filter(!is.na(logFC)) |>
  dplyr::select(
    logFC, cell, is_serine, sites_per_protein, pI_7mer, hydrophobicity_7mer,
    log_protein_abundance, pLDDT, pPSE_24, pPSE_12, is_IDR, secondary_structure_simple
  ) |>
  mutate(
    cell = factor(cell, levels = c("HEK293T", "HepG2", "Jurkat")),
    secondary_structure_simple = factor(secondary_structure_simple),
    is_serine = factor(is_serine, levels = c(0, 1), labels = c("Thr", "Ser")),
    is_IDR = factor(is_IDR, levels = c(0, 1), labels = c("Structured", "IDR"))
  ) |>
  filter(complete.cases(pick(everything())))

cat("Combined data: n =", nrow(reg_data_combined), "sites with complete data\n")
cat("  HEK293T:", sum(reg_data_combined$cell == "HEK293T"), "\n")
cat("  HepG2:", sum(reg_data_combined$cell == "HepG2"), "\n")
cat("  Jurkat:", sum(reg_data_combined$cell == "Jurkat"), "\n\n")

# Run main effects model
cat("--- Running Main Effects Model ---\n")
model_main <- lm(logFC ~ cell + is_serine + sites_per_protein + pI_7mer +
                   hydrophobicity_7mer + log_protein_abundance + pLDDT +
                   pPSE_24 + pPSE_12 + is_IDR + secondary_structure_simple,
                 data = reg_data_combined)

cat("R-squared:", round(summary(model_main)$r.squared, 4), "\n")
cat("Adj R-squared:", round(summary(model_main)$adj.r.squared, 4), "\n\n")

# Extract coefficients
coef_main <- as.data.frame(summary(model_main)$coefficients)
coef_main$term <- rownames(coef_main)
colnames(coef_main)[1:4] <- c("estimate", "std_error", "t_value", "p_value")
coef_main <- coef_main |>
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ ".",
      TRUE ~ ""
    )
  ) |>
  dplyr::select(term, estimate, std_error, t_value, p_value, significance)

# Calculate standardized coefficients
predictor_sds <- reg_data_combined |>
  summarize(
    sites_per_protein = sd(sites_per_protein, na.rm = TRUE),
    pI_7mer = sd(pI_7mer, na.rm = TRUE),
    hydrophobicity_7mer = sd(hydrophobicity_7mer, na.rm = TRUE),
    log_protein_abundance = sd(log_protein_abundance, na.rm = TRUE),
    pLDDT = sd(pLDDT, na.rm = TRUE),
    pPSE_24 = sd(pPSE_24, na.rm = TRUE),
    pPSE_12 = sd(pPSE_12, na.rm = TRUE)
  )

coef_main_std <- coef_main |>
  mutate(
    predictor_sd = case_when(
      term == "sites_per_protein" ~ predictor_sds$sites_per_protein,
      term == "pI_7mer" ~ predictor_sds$pI_7mer,
      term == "hydrophobicity_7mer" ~ predictor_sds$hydrophobicity_7mer,
      term == "log_protein_abundance" ~ predictor_sds$log_protein_abundance,
      term == "pLDDT" ~ predictor_sds$pLDDT,
      term == "pPSE_24" ~ predictor_sds$pPSE_24,
      term == "pPSE_12" ~ predictor_sds$pPSE_12,
      TRUE ~ 1
    ),
    std_estimate = estimate * predictor_sd,
    std_error_scaled = std_error * predictor_sd
  )

cat("--- Standardized Coefficients (Effect of 1 SD change) ---\n")
coef_main_std |>
  filter(!grepl("Intercept|secondary_structure", term)) |>
  dplyr::select(term, estimate, predictor_sd, std_estimate, p_value, significance) |>
  arrange(desc(abs(std_estimate))) |>
  print()

write_csv(
  coef_main_std,
  paste0(source_file_path, 'site_features/regression_coefficients_standardized.csv')
)


# ==============================================================================
# PART 3: FIGURE GENERATION
# ==============================================================================

# --- Standalone imports (for running PART 3 independently) ---
if (!exists("source_file_path")) {
  source('data_source.R')
}
library(tidyverse)
library(ggpubr)

# Create output directory
dir.create(paste0(figure_file_path, 'Figure6'), showWarnings = FALSE)

# Load pre-computed features if not already in environment
if (!exists("features_analysis")) {
  cat("Loading pre-computed site features...\n")
  features_analysis <- read_csv(
    paste0(source_file_path, 'site_features/OGlcNAc_site_features.csv'),
    show_col_types = FALSE
  )
  cat("Loaded", nrow(features_analysis), "sites\n")
}

# Load regression coefficients if not already in environment
if (!exists("coef_main_std")) {
  cat("Loading regression coefficients...\n")
  coef_main_std <- read_csv(
    paste0(source_file_path, 'site_features/regression_coefficients_standardized.csv'),
    show_col_types = FALSE
  )
}

# Load structural features (for pLDDT profiles)
structural_features <- read_csv(
  paste0(source_file_path, 'site_features/alphafold_structural_features.csv'),
  show_col_types = FALSE
)

# Prepare regression data if not already in environment
if (!exists("reg_data_combined")) {
  cat("Preparing regression data...\n")
  reg_data_combined <- features_analysis |>
    filter(!is.na(logFC)) |>
    dplyr::select(
      logFC, cell, is_serine, sites_per_protein, pI_7mer, hydrophobicity_7mer,
      log_protein_abundance, pLDDT, pPSE_24, pPSE_12, is_IDR, secondary_structure_simple
    ) |>
    mutate(
      cell = factor(cell, levels = c("HEK293T", "HepG2", "Jurkat")),
      secondary_structure_simple = factor(secondary_structure_simple),
      is_serine = factor(is_serine, levels = c(0, 1), labels = c("Thr", "Ser")),
      is_IDR = factor(is_IDR, levels = c(0, 1), labels = c("Structured", "IDR"))
    ) |>
    filter(complete.cases(pick(everything())))
  cat("Prepared", nrow(reg_data_combined), "sites for regression plots\n")
}
# --- End standalone imports ---

# ------------------------------------------------------------------------------
# Figure 6A: logFC Distribution by Cell Type
# ------------------------------------------------------------------------------

Figure6A <- features_analysis |>
  ggplot(aes(x = logFC, fill = cell)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = colors_cell) +
  labs(
    x = expression(log[2](Tuni/Ctrl)),
    y = "Density",
    fill = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, -5, 0),
    axis.title = element_text(size = 9, color = "black"),
    axis.text = element_text(size = 9, color = "black")
  )

ggsave(
  paste0(figure_file_path, 'Figure6/Figure6A.pdf'),
  Figure6A, width = 2, height = 1.5
)
cat("\nFigure 6A saved.\n")

# ------------------------------------------------------------------------------
# Figure 6B: Secondary Structure Distribution (Donut Plot)
# ------------------------------------------------------------------------------

ss_counts <- features_analysis |>
  filter(!is.na(secondary_structure_simple)) |>
  group_by(secondary_structure_simple) |>
  summarize(count = n(), .groups = "drop") |>
  mutate(
    category = case_when(
      secondary_structure_simple == "STRN" ~ "Strand",
      secondary_structure_simple == "BEND" ~ "Bend",
      secondary_structure_simple == "unstructured" ~ "Unstructured"
    ),
    category = factor(category, levels = c("Strand", "Bend", "Unstructured")),
    percentage = count / sum(count) * 100
  ) |>
  arrange(category) |>
  mutate(
    ymax = cumsum(percentage),
    ymin = lag(ymax, default = 0),
    label_pos = (ymin + ymax) / 2,
    label_text = paste0(round(percentage, 1), "%"),
    label_y_adjusted = case_when(
      category == "Strand" ~ label_pos - 2,
      category == "Bend" ~ label_pos + 3,
      TRUE ~ label_pos
    )
  )

ss_colors <- c(
  "Strand" = color_palette[2],
  "Bend" = color_palette[3],
  "Unstructured" = color_palette[6]
)

Figure6B <- ggplot(ss_counts, aes(ymax = ymax, ymin = ymin, xmax = 5, xmin = 2.8, fill = category)) +
  geom_rect(color = "white", linewidth = 0.3) +
  geom_text(
    data = ss_counts |> filter(percentage >= 10),
    aes(x = 3.9, y = label_pos, label = label_text),
    size = 3.2, color = "black"
  ) +
  geom_text(
    data = ss_counts |> filter(category == "Strand"),
    aes(x = 5.8, y = label_y_adjusted, label = paste0(category, " (", label_text, ")")),
    size = 2.8, color = "black", hjust = 0
  ) +
  geom_text(
    data = ss_counts |> filter(category == "Bend"),
    aes(x = 5.8, y = label_y_adjusted, label = paste0(category, " (", label_text, ")")),
    size = 2.8, color = "black", hjust = 0
  ) +
  geom_segment(
    data = ss_counts |> filter(category == "Strand"),
    aes(x = 5, xend = 5.6, y = label_pos, yend = label_y_adjusted),
    color = "grey50", linewidth = 0.3
  ) +
  geom_segment(
    data = ss_counts |> filter(category == "Bend"),
    aes(x = 5, xend = 5.6, y = label_pos, yend = label_y_adjusted),
    color = "grey50", linewidth = 0.3
  ) +
  annotate("text", x = 0, y = 0, label = "Secondary\nStructure",
           size = 3.2, fontface = "bold") +
  coord_polar(theta = "y") +
  xlim(c(0, 6.5)) +
  scale_fill_manual(values = ss_colors) +
  labs(fill = NULL) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 8)
  )

ggsave(
  paste0(figure_file_path, 'Figure6/Figure6B.pdf'),
  Figure6B, width = 2, height = 2
)
cat("Figure 6B saved.\n")

# ------------------------------------------------------------------------------
# Figure 6C: Standardized Regression Coefficients
# ------------------------------------------------------------------------------

coef_plot_std <- coef_main_std |>
  filter(!grepl("Intercept|secondary_structure", term)) |>
  mutate(
    term = case_when(
      term == "cellHepG2" ~ "HepG2 (vs HEK293T)",
      term == "cellJurkat" ~ "Jurkat (vs HEK293T)",
      term == "is_serineSer" ~ "Serine (vs Thr)",
      term == "sites_per_protein" ~ "Sites per Protein",
      term == "pI_7mer" ~ "pI (7-mer)",
      term == "hydrophobicity_7mer" ~ "Hydrophobicity",
      term == "log_protein_abundance" ~ "Protein Abundance",
      term == "pLDDT" ~ "pLDDT (Confidence)",
      term == "pPSE_24" ~ "pPSE_24 (Full Sphere)",
      term == "pPSE_12" ~ "pPSE_12 (Side Chain)",
      term == "is_IDRIDR" ~ "IDR (vs Structured)",
      TRUE ~ term
    ),
    is_significant = p_value < 0.05
  )

Figure6C <- coef_plot_std |>
  ggplot(aes(x = reorder(term, abs(std_estimate)), y = std_estimate, fill = is_significant)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = std_estimate - 1.96 * std_error_scaled,
                    ymax = std_estimate + 1.96 * std_error_scaled),
                width = 0.15, linewidth = 0.3, color = "grey30") +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = color_palette[1], "FALSE" = color_palette[6]),
                    labels = c("TRUE" = "p < 0.05", "FALSE" = "p >= 0.05"),
                    breaks = c("TRUE", "FALSE")) +
  labs(
    x = "",
    y = expression(paste("Standardized Coefficient (Effect of 1 SD on ", log[2], "FC)")),
    fill = "Significance"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 6, color = "black"),
    legend.text = element_text(size = 5, color = "black"),
    legend.key.size = unit(0.2, "cm"),
    axis.title = element_text(size = 8, color = "black"),
    axis.text = element_text(size = 8, color = "black")
  )

ggsave(
  paste0(figure_file_path, 'Figure6/Figure6C.pdf'),
  Figure6C, width = 3, height = 2
)
cat("Figure 6C saved.\n")

# ------------------------------------------------------------------------------
# Figure 6D: IDR vs Structured Effect Across Cell Types
# ------------------------------------------------------------------------------

interaction_data <- reg_data_combined |>
  group_by(cell, is_IDR) |>
  summarize(
    mean_logFC = mean(logFC, na.rm = TRUE),
    se_logFC = sd(logFC, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

Figure6D <- interaction_data |>
  ggplot(aes(x = cell, y = mean_logFC, fill = is_IDR)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean_logFC - se_logFC, ymax = mean_logFC + se_logFC),
                position = position_dodge(width = 0.7), width = 0.15,
                linewidth = 0.3, color = "grey30") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("Structured" = color_palette[4], "IDR" = color_palette[1])) +
  labs(
    x = "",
    y = expression(paste("Mean ", log[2], "FC (Tuni/Ctrl)")),
    fill = "Region Type"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 6, color = "black"),
    legend.text = element_text(size = 5, color = "black"),
    legend.key.size = unit(0.2, "cm"),
    axis.title = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 9, color = "black")
  )

ggsave(
  paste0(figure_file_path, 'Figure6/Figure6D.pdf'),
  Figure6D, width = 2, height = 2
)
cat("Figure 6D saved.\n")

# ------------------------------------------------------------------------------
# Figure 6E: HYOU1 (Q9Y4L1) - Structured Region Example
# ------------------------------------------------------------------------------
# HYOU1: Hypoxia up-regulated protein 1 (GRP170/ORP150)
# ER chaperone, HSP70 family member
# Length: 999 aa
#
# O-GlcNAc Site: T871
#   - Located in Domain D (SBDα): α-helical lid domain (residues 716-999)
#   - logFC: +2.47 (5.5-fold increase with tunicamycin)
#   - pLDDT: 94.12 (very high confidence - structured)
#   - adj.P.Val: 1.5e-05 (highly significant)
#
# Domain Architecture:
#   - Domain A (NBD): 1-430 - Nucleotide binding domain
#   - Domain B (SBDβ): 431-600 - Substrate binding β-sheet domain
#   - Domain C: 601-715 - Acidic unstructured loop
#   - Domain D (SBDα): 716-999 - α-helical lid domain (contains T871)

# pLDDT profile for HYOU1
Q9Y4L1_pLDDT <- structural_features |>
  filter(protein_id == "Q9Y4L1") |>
  dplyr::select(position, pLDDT)

# O-GlcNAc site data (HEK293T)
Q9Y4L1_site_data <- tibble(
  site_number = 871,
  logFC = 2.465,
  region = "Structured",
  pLDDT_site = 94.12
)

# Domain annotations for HYOU1
Q9Y4L1_domain_data <- tribble(
  ~domain, ~start_position, ~end_position,
  "NBD", 1, 430,
 "SBDβ", 431, 600,
  "Acidic loop", 601, 715,
 "SBDα", 716, 999
)

Q9Y4L1_rect_data <- tribble(
  ~start, ~end, ~top, ~bottom,
  0, 999, 1.5, 2.0
)

Q9Y4L1_site_labels <- Q9Y4L1_site_data |>
  mutate(label = paste0("T", site_number))

# Top panel: Domain architecture with O-GlcNAc site
Figure6E_top <- ggplot() +
  xlim(c(0, 999)) +
  ylim(c(-0.5, 3.5)) +
  # Protein bar
  geom_rect(
    data = Q9Y4L1_rect_data,
    aes(xmin = start, xmax = end, ymin = bottom, ymax = top),
    fill = "gray80", color = "black", linewidth = 0.3, alpha = 0.8
  ) +
  # Domain regions
  geom_rect(
    data = Q9Y4L1_domain_data,
    aes(xmin = start_position, xmax = end_position,
        ymin = Q9Y4L1_rect_data$bottom, ymax = Q9Y4L1_rect_data$top,
        fill = domain),
    alpha = 0.7
  ) +
  # Vertical line connecting site to protein bar
  geom_segment(
    data = Q9Y4L1_site_labels,
    aes(x = site_number, xend = site_number,
        y = Q9Y4L1_rect_data$bottom, yend = logFC),
    linetype = "dotted", color = "grey50", linewidth = 0.3
  ) +
  # Site marker at protein bar
  geom_point(
    data = Q9Y4L1_site_labels,
    aes(x = site_number, y = Q9Y4L1_rect_data$top),
    shape = 25, size = 1.5, fill = "black", color = "black"
  ) +
  # Site label
  geom_text(
    data = Q9Y4L1_site_labels,
    aes(x = site_number, y = Q9Y4L1_rect_data$top + 0.3, label = label),
    size = 2.5, color = "black", angle = 90, hjust = 0
  ) +
  # Horizontal line at y=0
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  # logFC point
  geom_point(
    data = Q9Y4L1_site_data,
    aes(x = site_number, y = logFC),
    color = color_palette[4], size = 2
  ) +
  scale_fill_manual(
    name = "Domain",
    values = c(
      "NBD" = color_palette[2],
     "SBDβ" = color_palette[3],
      "Acidic loop" = color_palette[6],
     "SBDα" = color_palette[4]
    )
  ) +
  labs(
    x = "",
    y = expression(log[2](Tuni/Ctrl)),
    title = "HYOU1 (Structured)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 9, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 7, color = "black"),
    legend.text = element_text(size = 6, color = "black"),
    legend.key.size = unit(0.2, "cm"),
    axis.title = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(5, 5, 0, 5)
  )

# Bottom panel: pLDDT profile
Figure6E_bottom <- ggplot(Q9Y4L1_pLDDT, aes(x = position, y = pLDDT)) +
  geom_line(color = "black", linewidth = 0.4) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red", linewidth = 0.4) +
  # Mark O-GlcNAc site
  geom_point(
    data = Q9Y4L1_site_data,
    aes(x = site_number, y = pLDDT_site),
    color = color_palette[4], size = 1.5
  ) +
  coord_cartesian(xlim = c(0, 999)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
  labs(x = "Amino acid position", y = "pLDDT") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 9, color = "black"),
    axis.text = element_text(size = 9, color = "black"),
    legend.position = "none",
    plot.margin = margin(0, 5, 5, 5)
  )

Figure6E <- ggarrange(
  Figure6E_top, Figure6E_bottom,
  ncol = 1, nrow = 2,
  heights = c(2, 1),
  align = "v",
  common.legend = FALSE
)

ggsave(
  paste0(figure_file_path, 'Figure6/Figure6E.pdf'),
  Figure6E, width = 3, height = 2
)
cat("Figure 6E (HYOU1) saved.\n")

# ------------------------------------------------------------------------------
# Figure 6E Additional: Structured Region Candidate Sites (PyMOL 3D Structures)
# ------------------------------------------------------------------------------
# These 4 candidate sites are O-GlcNAc sites significantly upregulated in
# structured regions (pLDDT >= 50) with annotated protein domains.
# All are novel sites (not in OGlcNAc Atlas).
#
# PyMOL scripts to generate 3D structure figures:
#   1. Figure6E_PRDX6_Y89_pymol.py
#      - P30041 Y89 (Peroxiredoxin-6) - HEK293T
#      - logFC: +1.50, pLDDT: 98.25
#      - Domain: Thioredoxin (5-169)
#
#   2. Figure6E_PRDX6_T95_pymol.py
#      - P30041 T95 (Peroxiredoxin-6) - HepG2
#      - logFC: +1.51, pLDDT: 97.62
#      - Domain: Thioredoxin (5-169)
#
#   3. Figure6E_DDX50_Y492_pymol.py
#      - Q9BQ39 Y492 (DEAD-box helicase 50) - HepG2
#      - logFC: +0.95, pLDDT: 94.94
#      - Domain: Helicase C-terminal (380-524)
#
#   4. Figure6E_PWP2_T23_pymol.py
#      - Q15269 T23 (Periodic tryptophan protein 2) - HepG2
#      - logFC: +0.45, pLDDT: 92.69
#      - Domain: WD repeat 1 (12-50)
#
# Run PyMOL scripts:
#   /opt/homebrew/bin/pymol -c -q Figure6E_PRDX6_Y89_pymol.py
#   /opt/homebrew/bin/pymol -c -q Figure6E_PRDX6_T95_pymol.py
#   /opt/homebrew/bin/pymol -c -q Figure6E_DDX50_Y492_pymol.py
#   /opt/homebrew/bin/pymol -c -q Figure6E_PWP2_T23_pymol.py
#
# Output directory: Figures/Figure6/Figure6E_candidates/
# EThcD spectra PDFs are also in the same directory.

# ------------------------------------------------------------------------------
# Figure 6F: HOXA13 (P31271) - IDR Region Example (3D Structure by PyMOL)
# ------------------------------------------------------------------------------
# HOXA13: Homeobox protein Hox-A13
# Transcription factor involved in limb development
# Length: 388 aa
#
# O-GlcNAc Site: S199
#   - Located in disordered N-terminal region
#   - logFC: -0.019 (essentially no change with tunicamycin)
#   - pLDDT: 41.62 (low confidence - IDR)
#   - adj.P.Val: 0.85 (not significant)
#
# 3D structure generated by Figure6F_HOXA13_pymol.py
# The PyMOL script visualizes:
#   - AlphaFold structure colored by pLDDT confidence
#   - O-GlcNAc site S199 highlighted as cyan sphere
#   - Transparent surface overlay showing secondary structure

cat("\n=== Figure 6 Complete ===\n")
cat("Generated panels:\n")
cat("  - Figure6A.pdf: logFC distribution by cell type\n")
cat("  - Figure6B.pdf: Secondary structure distribution\n")
cat("  - Figure6C.pdf: Standardized regression coefficients\n")
cat("  - Figure6D.pdf: IDR vs Structured effect\n")
cat("  - Figure6E.pdf: HYOU1 T871 (Structured, logFC=+2.47, pLDDT=94.1)\n")
cat("  - Figure6F: HOXA13 S199 3D structure (IDR, logFC=-0.02, pLDDT=41.6) - see PyMOL script\n")
