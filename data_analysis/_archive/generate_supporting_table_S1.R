# generate_supporting_table_S1.R
# Generate Supporting Table S1: Identification of O-GlcNAcylated proteins
#
# This script creates an Excel file with 3 sheets (HEK293T, HepG2, Jurkat)
# containing O-GlcNAcylated protein identification data with proper formatting.

library(tidyverse)
library(openxlsx)

# =============================================================================
# File paths
# =============================================================================

# Input: Filtered O-GlcNAc data (to get list of O-GlcNAc proteins)
oglcnac_files <- c(
  "HEK293T" = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/filtered/OGlcNAc_HEK293T.csv",
  "HepG2" = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/filtered/OGlcNAc_HepG2.csv",
  "Jurkat" = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/filtered/OGlcNAc_Jurkat.csv"
)

# Input: Protein.tsv files (source of protein-level data)
protein_files <- c(
  "HEK293T" = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_HEK293T/OGlyco/EThcD_OPair_TMT_Search/OGlyco/protein.tsv",
  "HepG2" = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_HepG2/OGlyco/EThcD_OPair_TMT_Search/OGlyco/protein.tsv",
  "Jurkat" = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_Jurkat/OGlyco/EThcD_OPair_TMT_Search/OGlyco/protein.tsv"
)

# Output file
output_file <- "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/supporting_tables/new_version/supporting_table_S1.xlsx"

# =============================================================================
# Define styles
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

# Data row style for text columns (1, 2, 8): Times New Roman, 11, left-aligned
data_style_text <- createStyle(
  fontName = "Times New Roman",
  fontSize = 11,
  halign = "left",
  valign = "center"
)

# Data row style for numeric columns (3-7): Times New Roman, 11, center-aligned
data_style_numeric <- createStyle(
  fontName = "Times New Roman",
  fontSize = 11,
  halign = "center",
  valign = "center"
)

# =============================================================================
# Process data for each cell type
# =============================================================================

process_cell_type <- function(cell_type, oglcnac_file, protein_file) {

  cat("Processing", cell_type, "...\n")

  # Read filtered O-GlcNAc data to get unique protein IDs
  oglcnac_data <- read_csv(oglcnac_file, show_col_types = FALSE)
  oglcnac_proteins <- unique(oglcnac_data$Protein.ID)
  cat("  O-GlcNAc proteins:", length(oglcnac_proteins), "\n")

  # Read protein.tsv
  protein_data <- read_tsv(protein_file, show_col_types = FALSE)
  cat("  Total proteins in protein.tsv:", nrow(protein_data), "\n")

  # Filter to O-GlcNAc proteins only
  protein_filtered <- protein_data %>%
    filter(`Protein ID` %in% oglcnac_proteins)
  cat("  Filtered O-GlcNAc proteins:", nrow(protein_filtered), "\n")

  # Select and rename columns to match old table format
  result <- protein_filtered %>%
    select(
      `UniProt Accession` = `Protein ID`,
      `Gene Symbol` = Gene,
      `Total Peptide` = `Total Peptides`,
      `Unique Peptide` = `Unique Peptides`,
      Length = Length,
      Coverage_pct = Coverage,
      `Protein Annotation` = `Protein Description`
    ) %>%
    # Calculate Coverage (number of amino acids) and Coverage Percentage
    # protein.tsv Coverage is already a percentage (e.g., 1.63 means 1.63%)
    # Coverage = Length × (Coverage_pct / 100), as integer
    # Coverage Percentage = Coverage_pct formatted as "X.XX%"
    mutate(
      Coverage = as.integer(round(Length * Coverage_pct / 100)),
      `Coverage Percentage` = paste0(sprintf("%.2f", Coverage_pct), "%")
    ) %>%
    # Reorder columns to match old table
    select(
      `UniProt Accession`,
      `Gene Symbol`,
      `Total Peptide`,
      `Unique Peptide`,
      Length,
      Coverage,
      `Coverage Percentage`,
      `Protein Annotation`
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

for (i in seq_along(cell_types)) {
  cell_type <- cell_types[i]
  sheet_label <- sheet_labels[i]
  sheet_name <- paste0("Sheet", i)

  # Process data
  data <- process_cell_type(cell_type, oglcnac_files[cell_type], protein_files[cell_type])

  # Add worksheet

addWorksheet(wb, sheet_name)

  # Title text
  title_text <- paste0("Table S1", sheet_label, ". Identification of O-GlcNAcylated proteins in ", cell_type, " cells")

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
  addStyle(wb, sheet_name, data_style_text, rows = 3:n_rows, cols = c(1, 2, 8), gridExpand = TRUE)

  # Numeric columns (3-7)
  addStyle(wb, sheet_name, data_style_numeric, rows = 3:n_rows, cols = 3:7, gridExpand = TRUE)

  # Set column widths (approximate from old table)
  setColWidths(wb, sheet_name, cols = 1, widths = 16.3)   # UniProt Accession
  setColWidths(wb, sheet_name, cols = 2, widths = 14.5)   # Gene Symbol
  setColWidths(wb, sheet_name, cols = 3, widths = 11.8)   # Total Peptide
  setColWidths(wb, sheet_name, cols = 4, widths = 13.6)   # Unique Peptide
  setColWidths(wb, sheet_name, cols = 5, widths = 6.6)    # Length
  setColWidths(wb, sheet_name, cols = 6, widths = 8.5)    # Coverage
  setColWidths(wb, sheet_name, cols = 7, widths = 18.3)   # Coverage Percentage
  setColWidths(wb, sheet_name, cols = 8, widths = 100)    # Protein Annotation

  cat("  Sheet", i, "created with", nrow(data), "proteins\n\n")
}

# =============================================================================
# Save workbook
# =============================================================================

saveWorkbook(wb, output_file, overwrite = TRUE)
cat("Supporting Table S1 saved to:\n", output_file, "\n")
