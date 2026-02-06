# Figure 6 (Filtered): O-GlcNAc Site Structural Feature Analysis
#
# This script regenerates Figure 6 A-D after removing Level1b PSMs with
# site probability < 0.75 (low confidence HCD-based localizations).
#
# Filtering criteria:
#   - Keep: Level1 (all)
#   - Keep: Level1b with probability >= 0.75
#   - Remove: Level1b with probability < 0.75

# ==============================================================================
# PART 1: DATA PREPARATION WITH FILTERING
# ==============================================================================

source('data_source.R')
source('data_source_DE.R')

library(Peptides)
library(ggpubr)
library(Biostrings)

# ------------------------------------------------------------------------------
# Section 1.1: Function to parse site probability
# ------------------------------------------------------------------------------

parse_site_probability <- function(prob_str) {
  prob <- str_extract(prob_str, "[0-9.]+(?=\\]$)")
  as.numeric(prob)
}

# ------------------------------------------------------------------------------
# Section 1.2: Load and Filter Site Data (Remove Low Prob Level1b)
# ------------------------------------------------------------------------------

cat("=== Loading and Filtering Site Data ===\n")
cat("Removing Level1b PSMs with site probability < 0.75\n\n")

# Function to load and filter site data
load_and_filter_site_data <- function(file_path, cell_type) {
  data <- read_csv(file_path, show_col_types = FALSE)

  # Parse site probability
  data <- data %>%
    mutate(Site_Prob = parse_site_probability(Site.Probabilities))

  # Count before filtering
  n_before <- nrow(data)
  n_level1 <- sum(data$Confidence.Level == "Level1", na.rm = TRUE)
  n_level1b <- sum(data$Confidence.Level == "Level1b", na.rm = TRUE)
  n_level1b_low_prob <- sum(data$Confidence.Level == "Level1b" & data$Site_Prob < 0.75, na.rm = TRUE)

  # Filter: Keep Level1, and Level1b with prob >= 0.75
  data_filtered <- data %>%
    filter(
      Confidence.Level == "Level1" |
      (Confidence.Level == "Level1b" & Site_Prob >= 0.75)
    )

  n_after <- nrow(data_filtered)

  cat(cell_type, ":\n")
  cat("  Before filtering:", n_before, "PSMs\n")
  cat("    Level1:", n_level1, "\n")
  cat("    Level1b:", n_level1b, "(", n_level1b_low_prob, "with prob < 0.75)\n")
  cat("  After filtering:", n_after, "PSMs\n")
  cat("  Removed:", n_before - n_after, "low probability Level1b PSMs\n\n")

  return(data_filtered)
}

# Load and filter each cell type
OGlcNAc_site_HEK293T_filtered <- load_and_filter_site_data(
  paste0(source_file_path, 'site/OGlcNAc_site_HEK293T.csv'), "HEK293T"
)
OGlcNAc_site_HepG2_filtered <- load_and_filter_site_data(
  paste0(source_file_path, 'site/OGlcNAc_site_HepG2.csv'), "HepG2"
)
OGlcNAc_site_Jurkat_filtered <- load_and_filter_site_data(
  paste0(source_file_path, 'site/OGlcNAc_site_Jurkat.csv'), "Jurkat"
)

# ------------------------------------------------------------------------------
# Section 1.3: Get Unique Sites After Filtering
# ------------------------------------------------------------------------------

# Get unique site_indices that remain after filtering
sites_kept_HEK293T <- unique(OGlcNAc_site_HEK293T_filtered$site_index)
sites_kept_HepG2 <- unique(OGlcNAc_site_HepG2_filtered$site_index)
sites_kept_Jurkat <- unique(OGlcNAc_site_Jurkat_filtered$site_index)

cat("=== Unique Sites After Filtering ===\n")
cat("HEK293T:", length(sites_kept_HEK293T), "sites\n")
cat("HepG2:", length(sites_kept_HepG2), "sites\n")
cat("Jurkat:", length(sites_kept_Jurkat), "sites\n\n")

# ------------------------------------------------------------------------------
# Section 1.4: Load Site-Level DE Data and Filter
# ------------------------------------------------------------------------------

OGlcNAc_site_DE_HEK293T <- read_csv(
  paste0(source_file_path, 'differential_analysis/OGlcNAc_site_DE_HEK293T.csv'),
  show_col_types = FALSE
) %>%
  filter(site_index %in% sites_kept_HEK293T)

OGlcNAc_site_DE_HepG2 <- read_csv(
  paste0(source_file_path, 'differential_analysis/OGlcNAc_site_DE_HepG2.csv'),
  show_col_types = FALSE
) %>%
  filter(site_index %in% sites_kept_HepG2)

OGlcNAc_site_DE_Jurkat <- read_csv(
  paste0(source_file_path, 'differential_analysis/OGlcNAc_site_DE_Jurkat.csv'),
  show_col_types = FALSE
) %>%
  filter(site_index %in% sites_kept_Jurkat)

cat("=== Site-Level DE Data After Filtering ===\n")
cat("HEK293T:", nrow(OGlcNAc_site_DE_HEK293T), "sites\n")
cat("HepG2:", nrow(OGlcNAc_site_DE_HepG2), "sites\n")
cat("Jurkat:", nrow(OGlcNAc_site_DE_Jurkat), "sites\n\n")

# ------------------------------------------------------------------------------
# Section 1.5: Extract Basic Features (same as original)
# ------------------------------------------------------------------------------

extract_basic_features <- function(site_data, de_data) {
  site_info <- site_data %>%
    group_by(site_index) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    dplyr::select(
      site_index, Protein.ID, Gene, Peptide,
      Protein.Start, Protein.End, modified_residue, site_number
    )

  features <- de_data %>%
    left_join(site_info, by = c("site_index", "Protein.ID", "Gene")) %>%
    mutate(is_serine = as.integer(modified_residue == "S"))

  return(features)
}

features_HEK293T <- extract_basic_features(OGlcNAc_site_HEK293T_filtered, OGlcNAc_site_DE_HEK293T)
features_HepG2 <- extract_basic_features(OGlcNAc_site_HepG2_filtered, OGlcNAc_site_DE_HepG2)
features_Jurkat <- extract_basic_features(OGlcNAc_site_Jurkat_filtered, OGlcNAc_site_DE_Jurkat)

# ------------------------------------------------------------------------------
# Section 1.6: Add sites_per_protein
# ------------------------------------------------------------------------------

add_sites_per_protein <- function(features_df) {
  sites_count <- features_df %>%
    group_by(Protein.ID) %>%
    summarize(sites_per_protein = n(), .groups = "drop")

  features_df %>%
    left_join(sites_count, by = "Protein.ID")
}

features_HEK293T <- add_sites_per_protein(features_HEK293T)
features_HepG2 <- add_sites_per_protein(features_HepG2)
features_Jurkat <- add_sites_per_protein(features_Jurkat)

# ------------------------------------------------------------------------------
# Section 1.7: Extract 7-mer Sequence and Calculate Sequence Properties
# ------------------------------------------------------------------------------

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
  features_df <- features_df %>%
    rowwise() %>%
    mutate(seq_7mer = extract_7mer_from_protein(Protein.ID, site_number, seq_lookup)) %>%
    ungroup()

  features_df <- features_df %>%
    rowwise() %>%
    mutate(seq_props = list(calculate_sequence_properties(seq_7mer))) %>%
    ungroup() %>%
    mutate(
      pI_7mer = sapply(seq_props, function(x) x$pI),
      hydrophobicity_7mer = sapply(seq_props, function(x) x$hydrophobicity)
    ) %>%
    dplyr::select(-seq_props)

  return(features_df)
}

cat("Extracting 7-mer sequences and calculating properties...\n")
features_HEK293T <- add_sequence_features(features_HEK293T, protein_seq_lookup)
features_HepG2 <- add_sequence_features(features_HepG2, protein_seq_lookup)
features_Jurkat <- add_sequence_features(features_Jurkat, protein_seq_lookup)

# ------------------------------------------------------------------------------
# Section 1.8: Add Protein Abundance
# ------------------------------------------------------------------------------

add_protein_abundance <- function(features_df, wp_data) {
  abundance_data <- wp_data[, c("UniProt_Accession", "AveExpr")]
  names(abundance_data)[2] <- "log_protein_abundance"
  features_df %>%
    left_join(abundance_data, by = c("Protein.ID" = "UniProt_Accession"))
}

features_HEK293T <- add_protein_abundance(features_HEK293T, WP_protein_DE_HEK293T)
features_HepG2 <- add_protein_abundance(features_HepG2, WP_protein_DE_HepG2)
features_Jurkat <- add_protein_abundance(features_Jurkat, WP_protein_DE_Jurkat)

# ------------------------------------------------------------------------------
# Section 1.9: Add Structural Features from Cached File
# ------------------------------------------------------------------------------

structural_cache_file <- paste0(source_file_path, 'site_features/alphafold_structural_features.csv')

if (file.exists(structural_cache_file)) {
  cat("Loading cached structural features...\n")
  structural_features <- read_csv(structural_cache_file, show_col_types = FALSE)

  merge_structural_features <- function(features_df, structural_df) {
    features_df %>%
      left_join(
        structural_df %>%
          dplyr::select(protein_id, position, pLDDT, pPSE_24, pPSE_12,
                        pPSE_24_smooth10, is_IDR, secondary_structure_simple),
        by = c("Protein.ID" = "protein_id", "site_number" = "position")
      )
  }

  features_HEK293T <- merge_structural_features(features_HEK293T, structural_features)
  features_HepG2 <- merge_structural_features(features_HepG2, structural_features)
  features_Jurkat <- merge_structural_features(features_Jurkat, structural_features)
} else {
  cat("WARNING: Structural features file not found. Skipping structural features.\n")
  features_HEK293T <- features_HEK293T %>%
    mutate(pPSE_24 = NA_real_, pPSE_12 = NA_real_, pPSE_24_smooth10 = NA_real_,
           is_IDR = NA_integer_, secondary_structure_simple = NA_character_, pLDDT = NA_real_)
  features_HepG2 <- features_HepG2 %>%
    mutate(pPSE_24 = NA_real_, pPSE_12 = NA_real_, pPSE_24_smooth10 = NA_real_,
           is_IDR = NA_integer_, secondary_structure_simple = NA_character_, pLDDT = NA_real_)
  features_Jurkat <- features_Jurkat %>%
    mutate(pPSE_24 = NA_real_, pPSE_12 = NA_real_, pPSE_24_smooth10 = NA_real_,
           is_IDR = NA_integer_, secondary_structure_simple = NA_character_, pLDDT = NA_real_)
}

# ------------------------------------------------------------------------------
# Section 1.10: Combine All Features
# ------------------------------------------------------------------------------

features_combined <- bind_rows(
  features_HEK293T %>% mutate(cell = "HEK293T"),
  features_HepG2 %>% mutate(cell = "HepG2"),
  features_Jurkat %>% mutate(cell = "Jurkat")
)

features_analysis <- features_combined %>%
  dplyr::select(
    site_index, Protein.ID, Gene, cell, site_number, Peptide, seq_7mer,
    logFC, adj.P.Val,
    is_serine, sites_per_protein, pI_7mer, hydrophobicity_7mer, log_protein_abundance,
    pLDDT, pPSE_24, pPSE_12, pPSE_24_smooth10, is_IDR, secondary_structure_simple
  )

cat("\n=== Final Analysis Dataframe Summary (Filtered) ===\n")
cat("Total observations:", nrow(features_analysis), "\n")
cat("Complete cases:", sum(complete.cases(features_analysis)), "\n")

# Save filtered features
write_csv(
  features_analysis,
  paste0(source_file_path, 'site_features/OGlcNAc_site_features_filtered.csv')
)
cat("Filtered features saved.\n")


# ==============================================================================
# PART 2: REGRESSION ANALYSIS
# ==============================================================================

cat("\n=== Combined Regression with Cell Type Factor (Filtered Data) ===\n\n")

reg_data_combined <- features_analysis %>%
  filter(!is.na(logFC)) %>%
  dplyr::select(
    logFC, cell, is_serine, sites_per_protein, pI_7mer, hydrophobicity_7mer,
    log_protein_abundance, pLDDT, pPSE_24, pPSE_12, is_IDR, secondary_structure_simple
  ) %>%
  mutate(
    cell = factor(cell, levels = c("HEK293T", "HepG2", "Jurkat")),
    secondary_structure_simple = factor(secondary_structure_simple),
    is_serine = factor(is_serine, levels = c(0, 1), labels = c("Thr", "Ser")),
    is_IDR = factor(is_IDR, levels = c(0, 1), labels = c("Structured", "IDR"))
  ) %>%
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
coef_main <- coef_main %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      p_value < 0.1 ~ ".",
      TRUE ~ ""
    )
  ) %>%
  dplyr::select(term, estimate, std_error, t_value, p_value, significance)

# Calculate standardized coefficients
predictor_sds <- reg_data_combined %>%
  summarize(
    sites_per_protein = sd(sites_per_protein, na.rm = TRUE),
    pI_7mer = sd(pI_7mer, na.rm = TRUE),
    hydrophobicity_7mer = sd(hydrophobicity_7mer, na.rm = TRUE),
    log_protein_abundance = sd(log_protein_abundance, na.rm = TRUE),
    pLDDT = sd(pLDDT, na.rm = TRUE),
    pPSE_24 = sd(pPSE_24, na.rm = TRUE),
    pPSE_12 = sd(pPSE_12, na.rm = TRUE)
  )

coef_main_std <- coef_main %>%
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
coef_main_std %>%
  filter(!grepl("Intercept|secondary_structure", term)) %>%
  dplyr::select(term, estimate, predictor_sd, std_estimate, p_value, significance) %>%
  arrange(desc(abs(std_estimate))) %>%
  print()

write_csv(
  coef_main_std,
  paste0(source_file_path, 'site_features/regression_coefficients_standardized_filtered.csv')
)


# ==============================================================================
# PART 3: FIGURE GENERATION (A-D)
# ==============================================================================

# Create output directory
dir.create(paste0(figure_file_path, 'Figure6_filtered'), showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Figure 6A: logFC Distribution by Cell Type
# ------------------------------------------------------------------------------

Figure6A <- features_analysis %>%
  ggplot(aes(x = logFC, fill = cell)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = colors_cell) +
  labs(
    x = expression(log[2](Tuni/Ctrl)),
    y = "Density",
    fill = NULL,
    title = "Filtered (removed Level1b prob < 0.75)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 8, hjust = 0.5),
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
  paste0(figure_file_path, 'Figure6_filtered/Figure6A_filtered.pdf'),
  Figure6A, width = 2.5, height = 1.8
)
cat("\nFigure 6A (filtered) saved.\n")

# ------------------------------------------------------------------------------
# Figure 6B: Secondary Structure Distribution (Donut Plot)
# ------------------------------------------------------------------------------

ss_counts <- features_analysis %>%
  filter(!is.na(secondary_structure_simple)) %>%
  group_by(secondary_structure_simple) %>%
  summarize(count = n(), .groups = "drop") %>%
  mutate(
    category = case_when(
      secondary_structure_simple == "STRN" ~ "Strand",
      secondary_structure_simple == "BEND" ~ "Bend",
      secondary_structure_simple == "unstructured" ~ "Unstructured"
    ),
    category = factor(category, levels = c("Strand", "Bend", "Unstructured")),
    percentage = count / sum(count) * 100
  ) %>%
  arrange(category) %>%
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
    data = ss_counts %>% filter(percentage >= 10),
    aes(x = 3.9, y = label_pos, label = label_text),
    size = 3.2, color = "black"
  ) +
  geom_text(
    data = ss_counts %>% filter(category == "Strand"),
    aes(x = 5.8, y = label_y_adjusted, label = paste0(category, " (", label_text, ")")),
    size = 2.8, color = "black", hjust = 0
  ) +
  geom_text(
    data = ss_counts %>% filter(category == "Bend"),
    aes(x = 5.8, y = label_y_adjusted, label = paste0(category, " (", label_text, ")")),
    size = 2.8, color = "black", hjust = 0
  ) +
  geom_segment(
    data = ss_counts %>% filter(category == "Strand"),
    aes(x = 5, xend = 5.6, y = label_pos, yend = label_y_adjusted),
    color = "grey50", linewidth = 0.3
  ) +
  geom_segment(
    data = ss_counts %>% filter(category == "Bend"),
    aes(x = 5, xend = 5.6, y = label_pos, yend = label_y_adjusted),
    color = "grey50", linewidth = 0.3
  ) +
  annotate("text", x = 0, y = 0, label = "Secondary\nStructure",
           size = 3.2, fontface = "bold") +
  coord_polar(theta = "y") +
  xlim(c(0, 6.5)) +
  scale_fill_manual(values = ss_colors) +
  labs(fill = NULL, title = "Filtered") +
  theme_void() +
  theme(
    plot.title = element_text(size = 8, hjust = 0.5),
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 8)
  )

ggsave(
  paste0(figure_file_path, 'Figure6_filtered/Figure6B_filtered.pdf'),
  Figure6B, width = 2, height = 2.2
)
cat("Figure 6B (filtered) saved.\n")

# ------------------------------------------------------------------------------
# Figure 6C: Standardized Regression Coefficients
# ------------------------------------------------------------------------------

coef_plot_std <- coef_main_std %>%
  filter(!grepl("Intercept|secondary_structure", term)) %>%
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

Figure6C <- coef_plot_std %>%
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
    fill = "Significance",
    title = "Filtered"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 8, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 6, color = "black"),
    legend.text = element_text(size = 5, color = "black"),
    legend.key.size = unit(0.2, "cm"),
    axis.title = element_text(size = 8, color = "black"),
    axis.text = element_text(size = 8, color = "black")
  )

ggsave(
  paste0(figure_file_path, 'Figure6_filtered/Figure6C_filtered.pdf'),
  Figure6C, width = 3.2, height = 2.2
)
cat("Figure 6C (filtered) saved.\n")

# ------------------------------------------------------------------------------
# Figure 6D: IDR vs Structured Effect Across Cell Types
# ------------------------------------------------------------------------------

interaction_data <- reg_data_combined %>%
  group_by(cell, is_IDR) %>%
  summarize(
    mean_logFC = mean(logFC, na.rm = TRUE),
    se_logFC = sd(logFC, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

Figure6D <- interaction_data %>%
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
    fill = "Region Type",
    title = "Filtered"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 8, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 6, color = "black"),
    legend.text = element_text(size = 5, color = "black"),
    legend.key.size = unit(0.2, "cm"),
    axis.title = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 9, color = "black")
  )

ggsave(
  paste0(figure_file_path, 'Figure6_filtered/Figure6D_filtered.pdf'),
  Figure6D, width = 2.2, height = 2.2
)
cat("Figure 6D (filtered) saved.\n")

# ==============================================================================
# Summary
# ==============================================================================

cat("\n=== Figure 6 (Filtered) Complete ===\n")
cat("Output directory:", paste0(figure_file_path, 'Figure6_filtered/'), "\n")
cat("Generated panels:\n")
cat("  - Figure6A_filtered.pdf: logFC distribution by cell type\n")
cat("  - Figure6B_filtered.pdf: Secondary structure distribution\n")
cat("  - Figure6C_filtered.pdf: Standardized regression coefficients\n")
cat("  - Figure6D_filtered.pdf: IDR vs Structured effect\n")
