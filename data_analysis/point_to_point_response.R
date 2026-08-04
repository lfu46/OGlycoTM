# Point-to-Point Response Script
# This script generates supplementary data for reviewer responses

# Import packages
library(tidyverse)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Source file path
source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'

# Output file path for point-to-point response
output_path <- paste0(source_file_path, 'point_to_point_response/')

# Create output directory if it doesn't exist
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
  cat("Created directory:", output_path, "\n")
}


# O-GlcNAc Level 1/1b PSM extraction with full search results ------------------
# This section extracts O-GlcNAc PSMs with confidence Level 1 or 1b from each
# cell type and merges with full PSM file columns for spectral annotation.
# Output files contain all search parameters needed for spectrum visualization.

cat("\n=== O-GlcNAc Level 1/1b PSM extraction for spectral annotation ===\n\n")

# Define O-GlcNAc compositions
OGlcNAc_compositions <- c(
  'HexNAt(1) % 299.1230',
  'HexNAt(1)TMT6plex(1) % 528.2859'
)

# Define PSM file paths for each cell type
psm_file_paths <- list(

HEK293T = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_HEK293T/OGlyco/EThcD_OPair_TMT_Search/OGlyco/psm.tsv',
  HepG2 = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_HepG2/OGlyco/EThcD_OPair_TMT_Search/OGlyco/psm.tsv',
  Jurkat = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_Jurkat/OGlyco/EThcD_OPair_TMT_Search/OGlyco/psm.tsv'
)

# Define bonafide file paths
bonafide_file_paths <- list(
  HEK293T = paste0(source_file_path, 'filtered/OGlyco_HEK293T_bonafide.csv'),
  HepG2 = paste0(source_file_path, 'filtered/OGlyco_HepG2_bonafide.csv'),
  Jurkat = paste0(source_file_path, 'filtered/OGlyco_Jurkat_bonafide.csv')
)

# Helper function to normalize column names for comparison
normalize_col_name <- function(x) gsub("[^a-zA-Z0-9]", "", tolower(x))

# Process each cell type
for (cell_type in names(psm_file_paths)) {
  cat("\n--- Processing", cell_type, "---\n")

  # Load bonafide data
  bonafide_df <- read_csv(bonafide_file_paths[[cell_type]], show_col_types = FALSE)
  cat("Loaded bonafide file:", nrow(bonafide_df), "rows,", ncol(bonafide_df), "columns\n")

  # Filter for O-GlcNAc Level 1/1b
  OGlcNAc_Level1 <- bonafide_df |>
    filter(Total.Glycan.Composition %in% OGlcNAc_compositions) |>
    filter(Confidence.Level %in% c("Level1", "Level1b"))

  cat("O-GlcNAc Level 1/1b:", nrow(OGlcNAc_Level1), "PSMs\n")
  cat("  - Level1:", sum(OGlcNAc_Level1$Confidence.Level == "Level1"), "\n")
  cat("  - Level1b:", sum(OGlcNAc_Level1$Confidence.Level == "Level1b"), "\n")

  # Load full PSM file
  cat("Loading PSM file...\n")
  psm_df <- read_tsv(psm_file_paths[[cell_type]], show_col_types = FALSE)
  cat("PSM file loaded:", nrow(psm_df), "rows,", ncol(psm_df), "columns\n")

  # Find PSM columns not already in bonafide (keep Spectrum for joining)
  bonafide_cols_normalized <- normalize_col_name(colnames(bonafide_df))
  psm_cols_normalized <- normalize_col_name(colnames(psm_df))
  new_cols_idx <- !psm_cols_normalized %in% bonafide_cols_normalized | colnames(psm_df) == "Spectrum"
  new_cols <- colnames(psm_df)[new_cols_idx]

  cat("Additional columns from PSM:", length(new_cols) - 1, "\n")

  # Select only new columns from PSM
  psm_subset <- psm_df |>
    select(all_of(new_cols))

  # Merge
  OGlcNAc_Level1_merged <- OGlcNAc_Level1 |>
    left_join(psm_subset, by = "Spectrum")

  cat("After merge:", nrow(OGlcNAc_Level1_merged), "rows,",
      ncol(OGlcNAc_Level1_merged), "columns\n")

  # Verify merge success
  n_matched <- sum(!is.na(OGlcNAc_Level1_merged$`Parent Scan Number`))
  cat("Successfully matched:", n_matched, "/", nrow(OGlcNAc_Level1_merged), "\n")

  # Save output
  output_file <- paste0(output_path, 'OGlcNAc_Level1_', cell_type, '_with_PSM_info.csv')
  write_csv(OGlcNAc_Level1_merged, file = output_file)
  cat("Saved to:", output_file, "\n")
}

cat("\n=== O-GlcNAc Level 1/1b PSM extraction complete ===\n")


# Composite Spectral Quality Score Calculation ---------------------------------
# This section calculates a composite score for ranking spectra by quality
# for selecting candidates for spectral annotation and visualization.
#
# Metrics used (with weights):
#   - O.Pair.Score (3): Primary glycan localization score
#   - Hyperscore (3): Core spectral match quality
#   - Delta.Score (2): Hyperscore - Nextscore, confidence of unique match
#   - Localization.Delta.Score (2): Glycan site localization confidence
#   - Probability (1): PeptideProphet confidence
#   - Expectation (1): E-value (inverted - lower is better)
#   - Purity (1): Spectrum cleanliness
#
# Each metric is converted to percentile rank (0-1), weighted, and summed.

cat("\n=== Calculating Composite Spectral Quality Score ===\n\n")

# Define input files (from previous step)
psm_info_files <- list(
  HEK293T = paste0(output_path, 'OGlcNAc_Level1_HEK293T_with_PSM_info.csv'),
  HepG2 = paste0(output_path, 'OGlcNAc_Level1_HepG2_with_PSM_info.csv'),
  Jurkat = paste0(output_path, 'OGlcNAc_Level1_Jurkat_with_PSM_info.csv')
)

# Define weights for composite score
score_weights <- list(
  O.Pair.Score = 3,
  Hyperscore = 3,
  Delta.Score = 2,
  Localization.Delta.Score = 2,
  Probability = 1,
  Expectation = 1,
  Purity = 1
)
total_weight <- sum(unlist(score_weights))

# Process each cell type
for (cell_type in names(psm_info_files)) {
  cat("\n--- Processing", cell_type, "---\n")

  # Load data
  df <- read_csv(psm_info_files[[cell_type]], show_col_types = FALSE)
  cat("Loaded:", nrow(df), "PSMs\n")

  # Calculate Delta Score (Hyperscore - Nextscore)
  df <- df |>
    mutate(Delta.Score = Hyperscore - Nextscore)

  # Calculate percentile ranks (0-1 scale) for each metric
  df <- df |>
    mutate(
      # Higher is better - direct percentile rank
      O.Pair.Score_rank = percent_rank(O.Pair.Score),
      Hyperscore_rank = percent_rank(Hyperscore),
      Delta.Score_rank = percent_rank(Delta.Score),
      Localization.Delta.Score_rank = percent_rank(Localization.Delta.Score),
      Probability_rank = percent_rank(Probability),
      Purity_rank = percent_rank(Purity),
      # Lower is better - invert the rank
      Expectation_rank = 1 - percent_rank(Expectation)
    )

  # Calculate composite score (weighted average of ranks)
  df <- df |>
    mutate(
      Composite_Score = (
        score_weights$O.Pair.Score * O.Pair.Score_rank +
        score_weights$Hyperscore * Hyperscore_rank +
        score_weights$Delta.Score * Delta.Score_rank +
        score_weights$Localization.Delta.Score * Localization.Delta.Score_rank +
        score_weights$Probability * Probability_rank +
        score_weights$Expectation * Expectation_rank +
        score_weights$Purity * Purity_rank
      ) / total_weight
    )

  # Rank by composite score (1 = best)
  df <- df |>
    arrange(desc(Composite_Score)) |>
    mutate(Quality_Rank = row_number())

  # Print summary
  cat(sprintf("Composite Score: Min=%.4f, Median=%.4f, Max=%.4f\n",
              min(df$Composite_Score), median(df$Composite_Score), max(df$Composite_Score)))

  # Show top 5
  cat("\nTop 5 spectra:\n")
  df |>
    select(Quality_Rank, Gene, Peptide, Composite_Score, O.Pair.Score, Hyperscore) |>
    head(5) |>
    print()

  # Save ranked output (temporary, will add site info next)
  ranked_output_file <- paste0(output_path, 'OGlcNAc_Level1_', cell_type, '_ranked.csv')
  write_csv(df, file = ranked_output_file)
  cat("\nSaved to:", ranked_output_file, "\n")
}

cat("\n=== Composite Score Calculation Complete ===\n")


# Add site information to ranked files -----------------------------------------
# Merge site_index, site_number, modified_residue, peptide_site, glycan_type
# from the existing OGlcNAc_site files

cat("\n=== Adding site information to ranked files ===\n\n")

# Define site file paths
site_files <- list(
  HEK293T = paste0(source_file_path, 'site/OGlcNAc_site_HEK293T.csv'),
  HepG2 = paste0(source_file_path, 'site/OGlcNAc_site_HepG2.csv'),
  Jurkat = paste0(source_file_path, 'site/OGlcNAc_site_Jurkat.csv')
)

# Site columns to add
site_cols_to_add <- c("Spectrum", "peptide_site", "modified_residue", "site_number", "site_index", "glycan_type")

for (cell_type in names(site_files)) {
  cat("Adding site info to", cell_type, "...\n")

  # Load ranked file
  ranked_file <- paste0(output_path, 'OGlcNAc_Level1_', cell_type, '_ranked.csv')
  ranked_df <- read_csv(ranked_file, show_col_types = FALSE)

  # Load site file and select only site columns
  site_df <- read_csv(site_files[[cell_type]], show_col_types = FALSE) |>
    select(all_of(site_cols_to_add))

  # Join site info to ranked file
  ranked_with_site <- ranked_df |>
    left_join(site_df, by = "Spectrum")

  # Reorder columns - put site info near the front
  col_order <- c(
    "Spectrum", "Peptide", "Gene", "Protein.ID",
    "site_index", "site_number", "modified_residue", "peptide_site", "glycan_type",
    "Quality_Rank", "Composite_Score"
  )
  other_cols <- setdiff(colnames(ranked_with_site), col_order)
  ranked_with_site <- ranked_with_site |>
    select(all_of(col_order), all_of(other_cols))

  # Save updated file
  write_csv(ranked_with_site, file = ranked_file)
  cat("  Added", sum(!is.na(ranked_with_site$site_index)), "site annotations\n")
}

cat("\n=== Site information added successfully ===\n")


# T cell GO terms (Figure 2A - Jurkat unique proteins) -------------------------
# This section extracts:
# 1. Unique O-GlcNAc proteins from Jurkat cells
# 2. Proteins from the 3 GO terms shown in Figure 2A for Jurkat:
#    - Leukocyte activation involved in immune response
#    - Leukocyte cell-cell adhesion
#    - T cell activation

cat("\n=== T cell GO terms (Figure 2A) ===\n\n")

# Load bonafide glyco data
OGlyco_HEK293T_bonafide <- read_csv(

  paste0(source_file_path, 'filtered/OGlyco_HEK293T_bonafide.csv'),
  show_col_types = FALSE
)
OGlyco_HepG2_bonafide <- read_csv(
  paste0(source_file_path, 'filtered/OGlyco_HepG2_bonafide.csv'),
  show_col_types = FALSE
)
OGlyco_Jurkat_bonafide <- read_csv(
  paste0(source_file_path, 'filtered/OGlyco_Jurkat_bonafide.csv'),
  show_col_types = FALSE
)

# Define O-GlcNAc compositions
OGlcNAc_compositions <- c(
  'HexNAt(1) % 299.1230',
  'HexNAt(1)TMT6plex(1) % 528.2859'
)

# Extract O-GlcNAc proteins from each cell type
OGlcNAc_proteins_HEK293T <- OGlyco_HEK293T_bonafide |>
  filter(Total.Glycan.Composition %in% OGlcNAc_compositions) |>
  pull(Protein.ID) |>
  unique()

OGlcNAc_proteins_HepG2 <- OGlyco_HepG2_bonafide |>
  filter(Total.Glycan.Composition %in% OGlcNAc_compositions) |>
  pull(Protein.ID) |>
  unique()

OGlcNAc_proteins_Jurkat <- OGlyco_Jurkat_bonafide |>
  filter(Total.Glycan.Composition %in% OGlcNAc_compositions) |>
  pull(Protein.ID) |>
  unique()

# Calculate unique Jurkat O-GlcNAc proteins
unique_OGlcNAc_Jurkat <- setdiff(
  setdiff(OGlcNAc_proteins_Jurkat, OGlcNAc_proteins_HEK293T),
  OGlcNAc_proteins_HepG2
)

cat("Total O-GlcNAc proteins in Jurkat:", length(OGlcNAc_proteins_Jurkat), "\n")
cat("Unique O-GlcNAc proteins in Jurkat:", length(unique_OGlcNAc_Jurkat), "\n\n")

# Convert UniProt IDs to gene symbols
convert_uniprot_to_symbol <- function(uniprot_ids) {
  symbols <- mapIds(
    org.Hs.eg.db,
    keys = uniprot_ids,
    column = "SYMBOL",
    keytype = "UNIPROT",
    multiVals = "first"
  )
  return(symbols)
}

# Get gene symbols for unique Jurkat proteins
unique_Jurkat_symbols <- convert_uniprot_to_symbol(unique_OGlcNAc_Jurkat)

# Create data frame with UniProt ID and Gene Symbol
unique_Jurkat_df <- data.frame(
  UniProt_ID = unique_OGlcNAc_Jurkat,
  Gene_Symbol = unique_Jurkat_symbols[unique_OGlcNAc_Jurkat],
  stringsAsFactors = FALSE
) |>
  arrange(Gene_Symbol)

# Export unique Jurkat proteins
write_csv(
  unique_Jurkat_df,
  file = paste0(output_path, 'Figure2A_unique_OGlcNAc_Jurkat_proteins.csv')
)
cat("Exported:", paste0(output_path, 'Figure2A_unique_OGlcNAc_Jurkat_proteins.csv'), "\n")

# --- GO Term Protein Lists ---

# Define the 3 GO terms from Figure 2A Jurkat panel
# Data extracted from enrichment/Figure2C_unique_Jurkat_OGlcNAc_GO.csv

GO_terms_Jurkat <- list(
  "Leukocyte_activation_immune_response" = list(
    GO_ID = "GO:0002366",
    Description = "leukocyte activation involved in immune response",
    pvalue = 8.169841172930072e-4,
    proteins = c("P15153", "P31146", "P08575", "P13796", "Q9HBG7",
                 "O15530", "P05107", "P02786", "Q9UJU2", "P17275", "Q96BY6")
  ),
  "Leukocyte_cell_cell_adhesion" = list(
    GO_ID = "GO:0007159",
    Description = "leukocyte cell-cell adhesion",
    pvalue = 0.001027866536388558,
    proteins = c("P15153", "Q13951", "P98172", "P31146", "P08575",
                 "O15530", "P05107", "Q01196", "Q9Y2R2", "Q92854",
                 "P02786", "P60709", "Q9UJU2")
  ),
  "T_cell_activation" = list(
    GO_ID = "GO:0042110",
    Description = "T cell activation",
    pvalue = 0.001355007230345075,
    proteins = c("P15153", "Q13951", "P98172", "P31146", "P08575",
                 "P13796", "Q9HBG7", "O15530", "Q01196", "Q9Y2R2",
                 "Q9C0K0", "P02786", "P60709", "Q9UJU2", "P17275")
  )
)

# Process each GO term and create output
all_GO_proteins <- data.frame()

for (term_name in names(GO_terms_Jurkat)) {
  term_data <- GO_terms_Jurkat[[term_name]]
  proteins <- term_data$proteins

  # Get gene symbols
  symbols <- convert_uniprot_to_symbol(proteins)

  # Create data frame
  term_df <- data.frame(
    GO_ID = term_data$GO_ID,
    GO_Description = term_data$Description,
    pvalue = term_data$pvalue,
    UniProt_ID = proteins,
    Gene_Symbol = symbols[proteins],
    stringsAsFactors = FALSE
  )

  all_GO_proteins <- bind_rows(all_GO_proteins, term_df)

  cat("\n", term_data$Description, "(", term_data$GO_ID, ")\n")
  cat("  p-value:", format(term_data$pvalue, scientific = TRUE), "\n")
  cat("  Count:", length(proteins), "\n")
  cat("  Genes:", paste(symbols[proteins], collapse = ", "), "\n")
}

# Export combined GO term protein list
write_csv(
  all_GO_proteins,
  file = paste0(output_path, 'Figure2A_Jurkat_GO_term_proteins.csv')
)
cat("\nExported:", paste0(output_path, 'Figure2A_Jurkat_GO_term_proteins.csv'), "\n")

# Create a summary table
GO_summary <- data.frame(
  GO_ID = sapply(GO_terms_Jurkat, function(x) x$GO_ID),
  Description = sapply(GO_terms_Jurkat, function(x) x$Description),
  pvalue = sapply(GO_terms_Jurkat, function(x) x$pvalue),
  Count = sapply(GO_terms_Jurkat, function(x) length(x$proteins)),
  Proteins = sapply(GO_terms_Jurkat, function(x) {
    symbols <- convert_uniprot_to_symbol(x$proteins)
    paste(symbols, collapse = ", ")
  }),
  stringsAsFactors = FALSE
)

write_csv(
  GO_summary,
  file = paste0(output_path, 'Figure2A_Jurkat_GO_term_summary.csv')
)
cat("Exported:", paste0(output_path, 'Figure2A_Jurkat_GO_term_summary.csv'), "\n")

cat("\n=== T cell GO terms section complete ===\n")


# ER O-GlcNAc proteins in HepG2 (Figure 5D) ------------------------------------
# This section extracts O-GlcNAc proteins that are:
# 1. Located in Endoplasmic Reticulum (ER)
# 2. From HepG2 cells
# 3. Significantly downregulated (logFC < 0 and adj.P.Val < 0.05)

cat("\n\n=== ER O-GlcNAc proteins in HepG2 (Figure 5D) ===\n\n")

# Load the combined location data with differential expression results
OGlcNAc_location_combined <- read_csv(
  paste0(source_file_path, 'subcellular_location/OGlcNAc_protein_location_combined.csv'),
  show_col_types = FALSE
)

cat("Loaded: OGlcNAc_protein_location_combined.csv\n")
cat("Total proteins with location annotation:", nrow(OGlcNAc_location_combined), "\n\n")

# Filter for ER proteins in HepG2
ER_HepG2_all <- OGlcNAc_location_combined |>
  filter(cell == "HepG2", Location == "Endoplasmic reticulum")

cat("All ER O-GlcNAc proteins in HepG2:", nrow(ER_HepG2_all), "\n")

# Filter for significantly downregulated proteins
# Criteria: logFC < 0 AND adj.P.Val < 0.05
ER_HepG2_sig_down <- ER_HepG2_all |>
  filter(logFC < 0, adj.P.Val < 0.05) |>
  arrange(logFC)  # Sort by logFC (most downregulated first)

cat("Significantly downregulated ER proteins (logFC < 0, adj.P.Val < 0.05):", nrow(ER_HepG2_sig_down), "\n\n")

# Display the results
if (nrow(ER_HepG2_sig_down) > 0) {
  cat("Significantly downregulated ER O-GlcNAc proteins in HepG2:\n")
  cat("-----------------------------------------------------------\n")

  # Print each protein
  for (i in 1:nrow(ER_HepG2_sig_down)) {
    row <- ER_HepG2_sig_down[i, ]
    cat(sprintf("%2d. %s (%s): logFC = %.3f, adj.P.Val = %.2e\n",
                i, row$Gene, row$Protein.ID, row$logFC, row$adj.P.Val))
  }

  # Export the list
  write_csv(
    ER_HepG2_sig_down,
    file = paste0(output_path, 'Figure5D_ER_HepG2_significantly_downregulated.csv')
  )
  cat("\nExported:", paste0(output_path, 'Figure5D_ER_HepG2_significantly_downregulated.csv'), "\n")

} else {
  cat("No significantly downregulated ER proteins found.\n")
  cat("Showing all ER proteins with negative logFC instead:\n\n")

  ER_HepG2_down <- ER_HepG2_all |>
    filter(logFC < 0) |>
    arrange(logFC)

  for (i in 1:nrow(ER_HepG2_down)) {
    row <- ER_HepG2_down[i, ]
    cat(sprintf("%2d. %s (%s): logFC = %.3f, adj.P.Val = %.2e\n",
                i, row$Gene, row$Protein.ID, row$logFC, row$adj.P.Val))
  }

  write_csv(
    ER_HepG2_down,
    file = paste0(output_path, 'Figure5D_ER_HepG2_downregulated_all.csv')
  )
  cat("\nExported:", paste0(output_path, 'Figure5D_ER_HepG2_downregulated_all.csv'), "\n")
}

# Also export all ER proteins in HepG2 for reference
write_csv(
  ER_HepG2_all |> arrange(logFC),
  file = paste0(output_path, 'Figure5D_ER_HepG2_all_proteins.csv')
)
cat("Exported:", paste0(output_path, 'Figure5D_ER_HepG2_all_proteins.csv'), "\n")

# Summary statistics
cat("\n--- Summary Statistics for ER O-GlcNAc proteins in HepG2 ---\n")
cat("Total ER proteins:", nrow(ER_HepG2_all), "\n")
cat("Downregulated (logFC < 0):", sum(ER_HepG2_all$logFC < 0), "\n")
cat("Upregulated (logFC > 0):", sum(ER_HepG2_all$logFC > 0), "\n")
cat("Significantly downregulated (logFC < 0, adj.P.Val < 0.05):", nrow(ER_HepG2_sig_down), "\n")
cat("Significantly upregulated (logFC > 0, adj.P.Val < 0.05):",
    sum(ER_HepG2_all$logFC > 0 & ER_HepG2_all$adj.P.Val < 0.05), "\n")
cat("Mean logFC:", round(mean(ER_HepG2_all$logFC), 3), "\n")
cat("Median logFC:", round(median(ER_HepG2_all$logFC), 3), "\n")

cat("\n=== ER O-GlcNAc proteins section complete ===\n")


# Mitochondria O-GlcNAc proteins in HepG2 (Figure 5E) --------------------------
# This section extracts O-GlcNAc proteins that are:
# 1. Located in Mitochondria
# 2. From HepG2 cells
# 3. Upregulated (logFC > 0) - with and without significance filter

cat("\n\n=== Mitochondria O-GlcNAc proteins in HepG2 (Figure 5E) ===\n\n")

# Filter for Mitochondria proteins in HepG2
Mito_HepG2_all <- OGlcNAc_location_combined |>
  filter(cell == "HepG2", Location == "Mitochondria")

cat("All Mitochondria O-GlcNAc proteins in HepG2:", nrow(Mito_HepG2_all), "\n")

# Filter for significantly upregulated proteins
# Criteria: logFC > 0 AND adj.P.Val < 0.05
Mito_HepG2_sig_up <- Mito_HepG2_all |>
  filter(logFC > 0, adj.P.Val < 0.05) |>
  arrange(desc(logFC))  # Sort by logFC (most upregulated first)

cat("Significantly upregulated Mito proteins (logFC > 0, adj.P.Val < 0.05):", nrow(Mito_HepG2_sig_up), "\n\n")

# Display significantly upregulated proteins
if (nrow(Mito_HepG2_sig_up) > 0) {
  cat("Significantly upregulated Mitochondria O-GlcNAc proteins in HepG2:\n")
  cat("-------------------------------------------------------------------\n")

  for (i in 1:nrow(Mito_HepG2_sig_up)) {
    row <- Mito_HepG2_sig_up[i, ]
    cat(sprintf("%2d. %s (%s): logFC = %.3f, adj.P.Val = %.2e\n",
                i, row$Gene, row$Protein.ID, row$logFC, row$adj.P.Val))
  }

  # Export the list
  write_csv(
    Mito_HepG2_sig_up,
    file = paste0(output_path, 'Figure5E_Mito_HepG2_significantly_upregulated.csv')
  )
  cat("\nExported:", paste0(output_path, 'Figure5E_Mito_HepG2_significantly_upregulated.csv'), "\n")
}

# Also show all upregulated (not just significant)
Mito_HepG2_up <- Mito_HepG2_all |>
  filter(logFC > 0) |>
  arrange(desc(logFC))

cat("\n\nAll upregulated Mitochondria O-GlcNAc proteins in HepG2 (logFC > 0):\n")
cat("---------------------------------------------------------------------\n")

for (i in 1:nrow(Mito_HepG2_up)) {
  row <- Mito_HepG2_up[i, ]
  sig_marker <- ifelse(row$adj.P.Val < 0.05, "*", "")
  cat(sprintf("%2d. %s (%s): logFC = %.3f, adj.P.Val = %.2e %s\n",
              i, row$Gene, row$Protein.ID, row$logFC, row$adj.P.Val, sig_marker))
}

write_csv(
  Mito_HepG2_up,
  file = paste0(output_path, 'Figure5E_Mito_HepG2_upregulated_all.csv')
)
cat("\nExported:", paste0(output_path, 'Figure5E_Mito_HepG2_upregulated_all.csv'), "\n")

# Export all Mito proteins in HepG2 for reference
write_csv(
  Mito_HepG2_all |> arrange(desc(logFC)),
  file = paste0(output_path, 'Figure5E_Mito_HepG2_all_proteins.csv')
)
cat("Exported:", paste0(output_path, 'Figure5E_Mito_HepG2_all_proteins.csv'), "\n")

# Summary statistics
cat("\n--- Summary Statistics for Mitochondria O-GlcNAc proteins in HepG2 ---\n")
cat("Total Mito proteins:", nrow(Mito_HepG2_all), "\n")
cat("Upregulated (logFC > 0):", sum(Mito_HepG2_all$logFC > 0), "\n")
cat("Downregulated (logFC < 0):", sum(Mito_HepG2_all$logFC < 0), "\n")
cat("Significantly upregulated (logFC > 0, adj.P.Val < 0.05):", nrow(Mito_HepG2_sig_up), "\n")
cat("Significantly downregulated (logFC < 0, adj.P.Val < 0.05):",
    sum(Mito_HepG2_all$logFC < 0 & Mito_HepG2_all$adj.P.Val < 0.05), "\n")
cat("Mean logFC:", round(mean(Mito_HepG2_all$logFC), 3), "\n")
cat("Median logFC:", round(median(Mito_HepG2_all$logFC), 3), "\n")

cat("\n=== Mitochondria O-GlcNAc proteins section complete ===\n")


# IDR O-GlcNAc sites with logFC slightly below 0 in HEK293T (Figure 6F) --------
# This section extracts O-GlcNAc sites that are:
# 1. In HEK293T cells
# 2. Located in IDR region (pLDDT < 50, is_IDR = 1)
# 3. logFC slightly lower than 0 (candidates for Figure 6F example)
#
# Current Figure 6F example: HOXA13 S199 (logFC = -0.019, pLDDT = 41.6)

cat("\n\n=== IDR O-GlcNAc sites with logFC ~ 0 in HEK293T (Figure 6F) ===\n\n")

# Load site features data
site_features <- read_csv(
  paste0(source_file_path, 'site_features/OGlcNAc_site_features.csv'),
  show_col_types = FALSE
)

cat("Loaded: OGlcNAc_site_features.csv\n")
cat("Total sites:", nrow(site_features), "\n\n")

# Filter for HEK293T cells
HEK293T_sites <- site_features |>
  filter(cell == "HEK293T")

cat("HEK293T sites:", nrow(HEK293T_sites), "\n")

# Filter for IDR regions (is_IDR = 1, which means pLDDT < 50)
HEK293T_IDR_sites <- HEK293T_sites |>
  filter(is_IDR == 1)

cat("HEK293T IDR sites (pLDDT < 50):", nrow(HEK293T_IDR_sites), "\n")

# Filter for logFC slightly below 0 (between -0.5 and 0)
# This captures sites with minimal change, slightly negative
HEK293T_IDR_slightly_negative <- HEK293T_IDR_sites |>
  filter(logFC < 0, logFC > -0.5) |>
  arrange(logFC) |>
  dplyr::select(
    site_index, Protein.ID, Gene, site_number,
    logFC, adj.P.Val, pLDDT, is_IDR,
    secondary_structure_simple, seq_7mer, Peptide
  )

cat("IDR sites with logFC between -0.5 and 0:", nrow(HEK293T_IDR_slightly_negative), "\n\n")

# Display the candidates
cat("Candidates for Figure 6F (IDR sites with logFC slightly < 0):\n")
cat("=============================================================\n")
cat("Current example: HOXA13 S199 (logFC = -0.019, pLDDT = 41.6)\n\n")

# Print all candidates sorted by logFC (closest to 0 first for easy comparison)
HEK293T_IDR_slightly_negative_display <- HEK293T_IDR_slightly_negative |>
  arrange(desc(logFC))  # Sort so closest to 0 is first

for (i in 1:min(nrow(HEK293T_IDR_slightly_negative_display), 50)) {
  row <- HEK293T_IDR_slightly_negative_display[i, ]
  sig_marker <- ifelse(row$adj.P.Val < 0.05, "*", "")
  cat(sprintf("%2d. %s %s%d: logFC = %.3f, pLDDT = %.1f, adj.P.Val = %.2e %s\n",
              i, row$Gene, substr(row$seq_7mer, 4, 4), row$site_number,
              row$logFC, row$pLDDT, row$adj.P.Val, sig_marker))
}

if (nrow(HEK293T_IDR_slightly_negative_display) > 50) {
  cat(sprintf("\n... and %d more sites (see CSV file for complete list)\n",
              nrow(HEK293T_IDR_slightly_negative_display) - 50))
}

# Export the list
write_csv(
  HEK293T_IDR_slightly_negative_display,
  file = paste0(output_path, 'Figure6F_HEK293T_IDR_sites_logFC_slightly_negative.csv')
)
cat("\nExported:", paste0(output_path, 'Figure6F_HEK293T_IDR_sites_logFC_slightly_negative.csv'), "\n")

# Also create a wider range list (logFC between -1 and 0) for more options
HEK293T_IDR_negative <- HEK293T_IDR_sites |>
  filter(logFC < 0, logFC > -1) |>
  arrange(desc(logFC)) |>
  dplyr::select(
    site_index, Protein.ID, Gene, site_number,
    logFC, adj.P.Val, pLDDT, is_IDR,
    secondary_structure_simple, seq_7mer, Peptide
  )

write_csv(
  HEK293T_IDR_negative,
  file = paste0(output_path, 'Figure6F_HEK293T_IDR_sites_logFC_negative_wider.csv')
)
cat("Exported:", paste0(output_path, 'Figure6F_HEK293T_IDR_sites_logFC_negative_wider.csv'),
    "(logFC between -1 and 0, n =", nrow(HEK293T_IDR_negative), ")\n")

# Summary statistics
cat("\n--- Summary Statistics for HEK293T IDR sites ---\n")
cat("Total IDR sites:", nrow(HEK293T_IDR_sites), "\n")
cat("IDR sites with logFC < 0:", sum(HEK293T_IDR_sites$logFC < 0, na.rm = TRUE), "\n")
cat("IDR sites with logFC > 0:", sum(HEK293T_IDR_sites$logFC > 0, na.rm = TRUE), "\n")
cat("IDR sites with -0.5 < logFC < 0:", nrow(HEK293T_IDR_slightly_negative), "\n")
cat("IDR sites with -1 < logFC < 0:", nrow(HEK293T_IDR_negative), "\n")
cat("Mean logFC (IDR sites):", round(mean(HEK293T_IDR_sites$logFC, na.rm = TRUE), 3), "\n")
cat("Median logFC (IDR sites):", round(median(HEK293T_IDR_sites$logFC, na.rm = TRUE), 3), "\n")

cat("\n=== IDR O-GlcNAc sites section complete ===\n")


# Figure 4C, 4D, 4E - Cell-type specific regulated O-GlcNAc proteins ------------
# This section extracts O-GlcNAc proteins for GO enrichment analysis:
# Figure 4C: Upregulated O-GlcNAc proteins in HEK293T (logFC > 0.5, adj.P.Val < 0.05)
# Figure 4D: Downregulated O-GlcNAc proteins in HepG2 (logFC < -0.5, adj.P.Val < 0.05)
# Figure 4E: Downregulated O-GlcNAc proteins in Jurkat (logFC < -0.5, adj.P.Val < 0.05)

# Load required packages
library(tidyverse)

# Define paths
source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'
output_path <- paste0(source_file_path, 'point_to_point_response/')

# Create output directory if it doesn't exist
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

cat("\n\n=== Cell-type specific regulated O-GlcNAc proteins (Figure 4C, 4D, 4E) ===\n\n")

# Load differential expression data
OGlcNAc_protein_DE_HEK293T <- read_csv(
  paste0(source_file_path, 'differential_analysis/OGlcNAc_protein_DE_HEK293T.csv'),
  show_col_types = FALSE
)
OGlcNAc_protein_DE_HepG2 <- read_csv(
  paste0(source_file_path, 'differential_analysis/OGlcNAc_protein_DE_HepG2.csv'),
  show_col_types = FALSE
)
OGlcNAc_protein_DE_Jurkat <- read_csv(
  paste0(source_file_path, 'differential_analysis/OGlcNAc_protein_DE_Jurkat.csv'),
  show_col_types = FALSE
)

cat("Loaded differential expression data for all cell types\n\n")

# --- Figure 4C: Upregulated O-GlcNAc proteins in HEK293T ---
cat("--- Figure 4C: Upregulated O-GlcNAc proteins in HEK293T ---\n")

OGlcNAc_up_HEK293T <- OGlcNAc_protein_DE_HEK293T |>
  filter(logFC > 0.5, adj.P.Val < 0.05) |>
  arrange(desc(logFC)) |>
  mutate(cell = "HEK293T", regulation = "Upregulated")

cat("Total upregulated proteins (logFC > 0.5, adj.P.Val < 0.05):", nrow(OGlcNAc_up_HEK293T), "\n\n")

# Display top proteins
cat("Top 20 upregulated O-GlcNAc proteins in HEK293T:\n")
for (i in 1:min(20, nrow(OGlcNAc_up_HEK293T))) {
  row <- OGlcNAc_up_HEK293T[i, ]
  cat(sprintf("%2d. %s (%s): logFC = %.3f, adj.P.Val = %.2e\n",
              i, row$Gene, row$Protein.ID, row$logFC, row$adj.P.Val))
}

# Export
write_csv(
  OGlcNAc_up_HEK293T,
  file = paste0(output_path, 'Figure4C_HEK293T_upregulated_OGlcNAc_proteins.csv')
)
cat("\nExported:", paste0(output_path, 'Figure4C_HEK293T_upregulated_OGlcNAc_proteins.csv'), "\n")

# --- Figure 4D: Downregulated O-GlcNAc proteins in HepG2 ---
cat("\n--- Figure 4D: Downregulated O-GlcNAc proteins in HepG2 ---\n")

OGlcNAc_down_HepG2 <- OGlcNAc_protein_DE_HepG2 |>
  filter(logFC < -0.5, adj.P.Val < 0.05) |>
  arrange(logFC) |>
  mutate(cell = "HepG2", regulation = "Downregulated")

cat("Total downregulated proteins (logFC < -0.5, adj.P.Val < 0.05):", nrow(OGlcNAc_down_HepG2), "\n\n")

# Display top proteins
cat("Top 20 downregulated O-GlcNAc proteins in HepG2:\n")
for (i in 1:min(20, nrow(OGlcNAc_down_HepG2))) {
  row <- OGlcNAc_down_HepG2[i, ]
  cat(sprintf("%2d. %s (%s): logFC = %.3f, adj.P.Val = %.2e\n",
              i, row$Gene, row$Protein.ID, row$logFC, row$adj.P.Val))
}

# Export
write_csv(
  OGlcNAc_down_HepG2,
  file = paste0(output_path, 'Figure4D_HepG2_downregulated_OGlcNAc_proteins.csv')
)
cat("\nExported:", paste0(output_path, 'Figure4D_HepG2_downregulated_OGlcNAc_proteins.csv'), "\n")

# --- Figure 4E: Downregulated O-GlcNAc proteins in Jurkat ---
cat("\n--- Figure 4E: Downregulated O-GlcNAc proteins in Jurkat ---\n")

OGlcNAc_down_Jurkat <- OGlcNAc_protein_DE_Jurkat |>
  filter(logFC < -0.5, adj.P.Val < 0.05) |>
  arrange(logFC) |>
  mutate(cell = "Jurkat", regulation = "Downregulated")

cat("Total downregulated proteins (logFC < -0.5, adj.P.Val < 0.05):", nrow(OGlcNAc_down_Jurkat), "\n\n")

# Display top proteins
cat("Top 20 downregulated O-GlcNAc proteins in Jurkat:\n")
for (i in 1:min(20, nrow(OGlcNAc_down_Jurkat))) {
  row <- OGlcNAc_down_Jurkat[i, ]
  cat(sprintf("%2d. %s (%s): logFC = %.3f, adj.P.Val = %.2e\n",
              i, row$Gene, row$Protein.ID, row$logFC, row$adj.P.Val))
}

# Export
write_csv(
  OGlcNAc_down_Jurkat,
  file = paste0(output_path, 'Figure4E_Jurkat_downregulated_OGlcNAc_proteins.csv')
)
cat("\nExported:", paste0(output_path, 'Figure4E_Jurkat_downregulated_OGlcNAc_proteins.csv'), "\n")

# Summary
cat("\n--- Summary ---\n")
cat("Figure 4C - HEK293T upregulated:", nrow(OGlcNAc_up_HEK293T), "proteins\n")
cat("Figure 4D - HepG2 downregulated:", nrow(OGlcNAc_down_HepG2), "proteins\n")
cat("Figure 4E - Jurkat downregulated:", nrow(OGlcNAc_down_Jurkat), "proteins\n")

cat("\n=== Cell-type specific regulated O-GlcNAc proteins section complete ===\n")
