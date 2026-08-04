# generate_supporting_table_S7.R
# Generate Supporting Table S7: Abundance changes of O-GlcNAcylation sites
#
# This script creates an Excel file with 3 sheets (HEK293T, HepG2, Jurkat)
# containing O-GlcNAcylation site differential expression data with proper formatting.
#
# Filtering criteria (site-level):
#   - Keep sites that have at least one high-confidence PSM
#   - High-confidence = Level1 (all) or Level1b with probability >= 0.75

library(tidyverse)
library(openxlsx)

# =============================================================================
# Helper function to parse site probability
# =============================================================================

parse_site_probability <- function(prob_str) {
  # Extract probability value from string like "S[0.95]" or "T[0.751]"
  prob <- str_extract(prob_str, "[0-9.]+(?=\\]$)")
  as.numeric(prob)
}

# =============================================================================
# File paths
# =============================================================================

source_file_path <- "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/"

# Input: O-GlcNAc site differential expression data
de_files <- c(
  "HEK293T" = paste0(source_file_path, "differential_analysis/OGlcNAc_site_DE_HEK293T.csv"),
  "HepG2" = paste0(source_file_path, "differential_analysis/OGlcNAc_site_DE_HepG2.csv"),
  "Jurkat" = paste0(source_file_path, "differential_analysis/OGlcNAc_site_DE_Jurkat.csv")
)

# Input: Site-level PSM data (for filtering by probability)
site_files <- c(
  "HEK293T" = paste0(source_file_path, "site/OGlcNAc_site_HEK293T.csv"),
  "HepG2" = paste0(source_file_path, "site/OGlcNAc_site_HepG2.csv"),
  "Jurkat" = paste0(source_file_path, "site/OGlcNAc_site_Jurkat.csv")
)

# Output file
output_file <- "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/supporting_tables/new_version/supporting_table_S7.xlsx"

# =============================================================================
# Define styles (same as S1-S6)
# =============================================================================

# Title row style: Times New Roman, 16, Bold, Blue (#0070C0), left-aligned
title_style <- createStyle(
  fontName = "Times New Roman",
  fontSize = 16,
  fontColour = "#0070C0",
  textDecoration = "bold",
  halign = "left",
  valign = "center"
)

# Header row style: Times New Roman, 12, Bold, Black, center-aligned
# With single black line on top, double black line on bottom
header_style <- createStyle(
  fontName = "Times New Roman",
  fontSize = 12,
  fontColour = "#000000",
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = c("top", "bottom"),
  borderStyle = c("thin", "double"),
  borderColour = c("#000000", "#000000")
)

# Data row style for text columns: Times New Roman, 11, left-aligned
data_style_text <- createStyle(
  fontName = "Times New Roman",
  fontSize = 11,
  halign = "left",
  valign = "center"
)

# Data row style for logFC column: Times New Roman, 11, center-aligned, 2 decimal places
data_style_logFC <- createStyle(
  fontName = "Times New Roman",
  fontSize = 11,
  halign = "center",
  valign = "center",
  numFmt = "0.00"
)

# Data row style for adj.P.Val column: Times New Roman, 11, center-aligned, scientific notation
data_style_pval <- createStyle(
  fontName = "Times New Roman",
  fontSize = 11,
  halign = "center",
  valign = "center",
  numFmt = "0.0E+00"
)

# =============================================================================
# Process data for each cell type
# =============================================================================

process_cell_type <- function(cell_type, de_file, site_file) {

  cat("Processing", cell_type, "...\n")

  # Read site data to get probability information
  site_data <- read_csv(site_file, show_col_types = FALSE)

  # Parse site probability and identify high-confidence sites
  site_data <- site_data %>%
    mutate(Site_Prob = parse_site_probability(Site.Probabilities))

  # Count before filtering
  n_level1 <- sum(site_data$Confidence.Level == "Level1", na.rm = TRUE)
  n_level1b <- sum(site_data$Confidence.Level == "Level1b", na.rm = TRUE)
  n_level1b_low_prob <- sum(site_data$Confidence.Level == "Level1b" & site_data$Site_Prob < 0.75, na.rm = TRUE)

  cat("  PSMs: Level1 =", n_level1, ", Level1b =", n_level1b, "(", n_level1b_low_prob, "with prob < 0.75)\n")

  # Get unique site_indices with at least one high-confidence PSM
  high_conf_sites <- site_data %>%
    filter(
      Confidence.Level == "Level1" |
      (Confidence.Level == "Level1b" & Site_Prob >= 0.75)
    ) %>%
    pull(site_index) %>%
    unique()

  # Read differential expression data
  de_data <- read_csv(de_file, show_col_types = FALSE)
  n_before <- nrow(de_data)

  # Filter to keep only high-confidence sites
  de_data_filtered <- de_data %>%
    filter(site_index %in% high_conf_sites)

  n_after <- nrow(de_data_filtered)

  cat("  Sites before filtering:", n_before, "\n")
  cat("  Sites after filtering:", n_after, "\n")
  cat("  Removed:", n_before - n_after, "sites\n")

  # Select and rename columns to match table format
  # Extract Site from site_index (e.g., P35658_T1845 -> T1845)
  result <- de_data_filtered %>%
    mutate(
      Site = sub("^[^_]+_", "", site_index)  # Remove everything before and including first underscore
    ) %>%
    select(
      `UniProt Accession` = Protein.ID,
      `Gene Symbol` = Gene,
      Site = Site,
      `Avg(log2(Tuni/Ctrl))` = logFC,
      `Adjusted P value` = adj.P.Val,
      `Protein Annotation` = Protein.Description
    ) %>%
    # Sort by UniProt Accession, then by Site
    arrange(`UniProt Accession`, Site)

  return(result)
}

# =============================================================================
# Create workbook
# =============================================================================

wb <- createWorkbook()

cell_types <- c("HEK293T", "HepG2", "Jurkat")
sheet_labels <- c("A", "B", "C")

# Define column types for formatting
# Text columns: 1 (UniProt Accession), 2 (Gene Symbol), 3 (Site), 6 (Protein Annotation)
text_cols <- c(1, 2, 3, 6)

for (i in seq_along(cell_types)) {
  cell_type <- cell_types[i]
  sheet_label <- sheet_labels[i]
  sheet_name <- paste0("Sheet", i)

  # Process data
  data <- process_cell_type(cell_type, de_files[cell_type], site_files[cell_type])

  # Add worksheet
  addWorksheet(wb, sheet_name)

  # Title text
  title_text <- paste0("Table S7", sheet_label, ". Abundance changes of O-GlcNAcylation sites upon tunicamycin treatment in ", cell_type, " cells")

  # Write title row (row 1)
  writeData(wb, sheet_name, title_text, startRow = 1, startCol = 1, colNames = FALSE)

  # Merge title row across all columns
  mergeCells(wb, sheet_name, cols = 1:6, rows = 1)

  # Write header row (row 2)
  writeData(wb, sheet_name, as.data.frame(t(colnames(data))), startRow = 2, startCol = 1, colNames = FALSE)

  # Write data rows (starting row 3)
  writeData(wb, sheet_name, data, startRow = 3, startCol = 1, colNames = FALSE)

  # Apply styles
  # Title style (row 1)
  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1, gridExpand = TRUE)

  # Header style (row 2, all columns)
  addStyle(wb, sheet_name, header_style, rows = 2, cols = 1:6, gridExpand = TRUE)

  # Data styles (row 3 onwards)
  n_rows <- nrow(data) + 2  # +2 for title and header rows

  # Text columns (1, 2, 3, 6)
  addStyle(wb, sheet_name, data_style_text, rows = 3:n_rows, cols = text_cols, gridExpand = TRUE)

  # logFC column (4) - 2 decimal places
  addStyle(wb, sheet_name, data_style_logFC, rows = 3:n_rows, cols = 4, gridExpand = TRUE)

  # adj.P.Val column (5) - scientific notation
  addStyle(wb, sheet_name, data_style_pval, rows = 3:n_rows, cols = 5, gridExpand = TRUE)

  # Set column widths
  setColWidths(wb, sheet_name, cols = 1, widths = 16)    # UniProt Accession
  setColWidths(wb, sheet_name, cols = 2, widths = 14)    # Gene Symbol
  setColWidths(wb, sheet_name, cols = 3, widths = 10)    # Site
  setColWidths(wb, sheet_name, cols = 4, widths = 20)    # Avg(log2(Tuni/Ctrl))
  setColWidths(wb, sheet_name, cols = 5, widths = 16)    # Adjusted P value
  setColWidths(wb, sheet_name, cols = 6, widths = 100)   # Protein Annotation

  cat("  Sheet", i, "created with", nrow(data), "sites\n\n")
}

# =============================================================================
# Save workbook
# =============================================================================

saveWorkbook(wb, output_file, overwrite = TRUE)
cat("Supporting Table S7 saved to:\n", output_file, "\n")
