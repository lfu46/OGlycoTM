# generate_supporting_table_S4.R
# Generate Supporting Table S4: Identification of total proteins
#
# This script creates an Excel file with 3 sheets (HEK293T, HepG2, Jurkat)
# containing whole proteome protein identification data with proper formatting.

library(tidyverse)
library(openxlsx)
library(Biostrings)  # For reading FASTA files

# =============================================================================
# File paths
# =============================================================================

# Input: WP filtered data (PSM-level)
wp_files <- c(
  "HEK293T" = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/filtered/WP_HEK293T_filtered.csv",
  "HepG2" = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/filtered/WP_HepG2_filtered.csv",
  "Jurkat" = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/filtered/WP_Jurkat_filtered.csv"
)

# Input: FASTA file for protein lengths
fasta_file <- "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/reference/uniprotkb_reviewed_true_AND_model_organ_2026_01_09.fasta"

# Output file
output_file <- "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/supporting_tables/new_version/supporting_table_S4.xlsx"

# =============================================================================
# Define styles (same as S1-S3)
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

# Data row style for numeric columns: Times New Roman, 11, center-aligned
data_style_numeric <- createStyle(
  fontName = "Times New Roman",
  fontSize = 11,
  halign = "center",
  valign = "center"
)

# Data row style for Coverage Percentage: Times New Roman, 11, center-aligned, 2 decimal places with %
data_style_percent <- createStyle(
  fontName = "Times New Roman",
  fontSize = 11,
  halign = "center",
  valign = "center",
  numFmt = "0.00%"
)

# =============================================================================
# Read FASTA file and extract protein lengths
# =============================================================================

cat("Reading FASTA file...\n")
fasta_seqs <- readAAStringSet(fasta_file)

# Extract UniProt accession from FASTA headers (format: >sp|P12345|NAME_HUMAN ...)
fasta_headers <- names(fasta_seqs)
uniprot_ids <- sub("^[a-z]+\\|([^|]+)\\|.*", "\\1", fasta_headers)

# Create protein length lookup table
protein_lengths <- data.frame(
  UniProt_Accession = uniprot_ids,
  Protein_Length = width(fasta_seqs),
  stringsAsFactors = FALSE
)
cat("  Loaded", nrow(protein_lengths), "protein sequences from FASTA\n\n")

# =============================================================================
# Helper function to extract peptide sequence from format "K.PEPTIDE.R"
# =============================================================================

extract_peptide_seq <- function(peptide) {
  # Remove flanking residues (everything before first . and after last .)
  seq <- sub("^[A-Z-]\\.", "", peptide)  # Remove prefix like "K."
  seq <- sub("\\.[A-Z-]$", "", seq)       # Remove suffix like ".R"
  # Remove modification symbols like * for length calculation
  seq <- gsub("[^A-Z]", "", seq)
  return(seq)
}

# =============================================================================
# Process data for each cell type
# =============================================================================

process_cell_type <- function(cell_type, wp_file, protein_lengths) {

  cat("Processing", cell_type, "...\n")

  # Read WP filtered data
  wp_data <- read_csv(wp_file, show_col_types = FALSE)
  cat("  Total PSMs:", nrow(wp_data), "\n")

  # Extract peptide sequences and calculate lengths
  wp_data <- wp_data %>%
    mutate(
      Peptide_Seq = sapply(Peptide, extract_peptide_seq),
      Peptide_Length = nchar(Peptide_Seq)
    )

  # Aggregate to protein level
  protein_data <- wp_data %>%
    group_by(UniProt_Accession, Gene.Symbol, Annotation) %>%
    summarise(
      Total_Peptide = n(),
      Unique_Peptide = n_distinct(Peptide_Seq),
      .groups = "drop"
    )

  # Calculate coverage from unique peptides
  coverage_data <- wp_data %>%
    select(UniProt_Accession, Peptide_Seq, Peptide_Length) %>%
    distinct() %>%
    group_by(UniProt_Accession) %>%
    summarise(Coverage = sum(Peptide_Length), .groups = "drop")

  # Join coverage data
  protein_data <- protein_data %>%
    left_join(coverage_data, by = "UniProt_Accession")

  cat("  Unique proteins:", nrow(protein_data), "\n")

  # Join with protein lengths from FASTA
  protein_data <- protein_data %>%
    left_join(protein_lengths, by = c("UniProt_Accession" = "UniProt_Accession"))

  # Calculate Coverage Percentage (as decimal for Excel formatting)
  protein_data <- protein_data %>%
    mutate(
      Coverage_Percentage = ifelse(!is.na(Protein_Length) & Protein_Length > 0,
                                    Coverage / Protein_Length,
                                    NA)
    )

  # Check for missing protein lengths
  missing_lengths <- sum(is.na(protein_data$Protein_Length))
  if (missing_lengths > 0) {
    cat("  Warning:", missing_lengths, "proteins missing length from FASTA\n")
  }

  # Select and rename columns to match table format
  result <- protein_data %>%
    select(
      `UniProt Accession` = UniProt_Accession,
      `Gene Symbol` = Gene.Symbol,
      `Total Peptide` = Total_Peptide,
      `Unique Peptide` = Unique_Peptide,
      Length = Protein_Length,
      Coverage = Coverage,
      `Coverage Percentage` = Coverage_Percentage,
      `Protein Annotation` = Annotation
    ) %>%
    # Sort by UniProt Accession
    arrange(`UniProt Accession`)

  return(result)
}

# =============================================================================
# Create workbook
# =============================================================================

wb <- createWorkbook()

cell_types <- c("HEK293T", "HepG2", "Jurkat")
sheet_labels <- c("A", "B", "C")

# Define column types for formatting
# Text columns: 1 (UniProt Accession), 2 (Gene Symbol), 8 (Protein Annotation)
# Numeric columns: 3 (Total Peptide), 4 (Unique Peptide), 5 (Length), 6 (Coverage)
# Percent column: 7 (Coverage Percentage)
text_cols <- c(1, 2, 8)
numeric_cols <- c(3, 4, 5, 6)

for (i in seq_along(cell_types)) {
  cell_type <- cell_types[i]
  sheet_label <- sheet_labels[i]
  sheet_name <- paste0("Sheet", i)

  # Process data
  data <- process_cell_type(cell_type, wp_files[cell_type], protein_lengths)

  # Add worksheet
  addWorksheet(wb, sheet_name)

  # Title text
  title_text <- paste0("Table S4", sheet_label, ". Identification of total proteins in ", cell_type, " cells")

  # Write title row (row 1)
  writeData(wb, sheet_name, title_text, startRow = 1, startCol = 1, colNames = FALSE)

  # Merge title row across all columns
  mergeCells(wb, sheet_name, cols = 1:8, rows = 1)

  # Write header row (row 2)
  writeData(wb, sheet_name, as.data.frame(t(colnames(data))), startRow = 2, startCol = 1, colNames = FALSE)

  # Write data rows (starting row 3)
  writeData(wb, sheet_name, data, startRow = 3, startCol = 1, colNames = FALSE)

  # Apply styles
  # Title style (row 1)
  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1, gridExpand = TRUE)

  # Header style (row 2, all columns)
  addStyle(wb, sheet_name, header_style, rows = 2, cols = 1:8, gridExpand = TRUE)

  # Data styles (row 3 onwards)
  n_rows <- nrow(data) + 2  # +2 for title and header rows

  # Text columns (1, 2, 8)
  addStyle(wb, sheet_name, data_style_text, rows = 3:n_rows, cols = text_cols, gridExpand = TRUE)

  # Numeric columns (3, 4, 5, 6)
  addStyle(wb, sheet_name, data_style_numeric, rows = 3:n_rows, cols = numeric_cols, gridExpand = TRUE)

  # Coverage Percentage column (7) - percentage format
  addStyle(wb, sheet_name, data_style_percent, rows = 3:n_rows, cols = 7, gridExpand = TRUE)

  # Set column widths
  setColWidths(wb, sheet_name, cols = 1, widths = 16)    # UniProt Accession
  setColWidths(wb, sheet_name, cols = 2, widths = 14)    # Gene Symbol
  setColWidths(wb, sheet_name, cols = 3, widths = 12)    # Total Peptide
  setColWidths(wb, sheet_name, cols = 4, widths = 14)    # Unique Peptide
  setColWidths(wb, sheet_name, cols = 5, widths = 8)     # Length
  setColWidths(wb, sheet_name, cols = 6, widths = 10)    # Coverage
  setColWidths(wb, sheet_name, cols = 7, widths = 18)    # Coverage Percentage
  setColWidths(wb, sheet_name, cols = 8, widths = 100)   # Protein Annotation

  cat("  Sheet", i, "created with", nrow(data), "proteins\n\n")
}

# =============================================================================
# Save workbook
# =============================================================================

saveWorkbook(wb, output_file, overwrite = TRUE)
cat("Supporting Table S4 saved to:\n", output_file, "\n")
