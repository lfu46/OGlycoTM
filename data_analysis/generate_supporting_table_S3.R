# generate_supporting_table_S3.R
# Generate Supporting Table S3: Abundance changes of O-GlcNAcylated proteins
#
# This script creates an Excel file with 3 sheets (HEK293T, HepG2, Jurkat)
# containing O-GlcNAcylated protein differential expression data with proper formatting.

library(tidyverse)
library(openxlsx)

# =============================================================================
# File paths
# =============================================================================

# Input: O-GlcNAc protein differential expression data
de_files <- c(
  "HEK293T" = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/differential_analysis/OGlcNAc_protein_DE_HEK293T.csv",
  "HepG2" = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/differential_analysis/OGlcNAc_protein_DE_HepG2.csv",
  "Jurkat" = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/differential_analysis/OGlcNAc_protein_DE_Jurkat.csv"
)

# Output file
output_file <- "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/supporting_tables/new_version/supporting_table_S3.xlsx"

# =============================================================================
# Define styles (same as S1, S2)
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

# Data row style for text columns (1, 2, 5): Times New Roman, 11, left-aligned
data_style_text <- createStyle(
  fontName = "Times New Roman",
  fontSize = 11,
  halign = "left",
  valign = "center"
)

# Data row style for logFC column (3): Times New Roman, 11, center-aligned, 2 decimal places
data_style_logFC <- createStyle(
  fontName = "Times New Roman",
  fontSize = 11,
  halign = "center",
  valign = "center",
  numFmt = "0.00"
)

# Data row style for adj.P.Val column (4): Times New Roman, 11, center-aligned, scientific notation
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

process_cell_type <- function(cell_type, de_file) {

  cat("Processing", cell_type, "...\n")

  # Read differential expression data
  de_data <- read_csv(de_file, show_col_types = FALSE)
  cat("  Total proteins:", nrow(de_data), "\n")

  # Select and rename columns to match table format
  # Keep numeric values as numbers (formatting applied via Excel styles)
  result <- de_data %>%
    select(
      `UniProt Accession` = Protein.ID,
      `Gene Symbol` = Gene,
      `Avg(log2(Tuni/Ctrl))` = logFC,
      `Adjusted P value` = adj.P.Val,
      `Protein Annotation` = Protein.Description
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
# Text columns: 1 (UniProt Accession), 2 (Gene Symbol), 5 (Protein Annotation)
# Numeric columns: 3 (logFC - 2 decimals), 4 (adj.P.Val - scientific notation)
text_cols <- c(1, 2, 5)

for (i in seq_along(cell_types)) {
  cell_type <- cell_types[i]
  sheet_label <- sheet_labels[i]
  sheet_name <- paste0("Sheet", i)

  # Process data
  data <- process_cell_type(cell_type, de_files[cell_type])

  # Add worksheet
  addWorksheet(wb, sheet_name)

  # Title text
  title_text <- paste0("Table S3", sheet_label, ". Abundance changes of O-GlcNAcylated proteins upon tunicamycin treatment in ", cell_type, " cells")

  # Write title row (row 1)
  writeData(wb, sheet_name, title_text, startRow = 1, startCol = 1, colNames = FALSE)

  # Merge title row across all columns
  mergeCells(wb, sheet_name, cols = 1:5, rows = 1)

  # Write header row (row 2)
  writeData(wb, sheet_name, as.data.frame(t(colnames(data))), startRow = 2, startCol = 1, colNames = FALSE)

  # Write data rows (starting row 3)
  writeData(wb, sheet_name, data, startRow = 3, startCol = 1, colNames = FALSE)

  # Apply styles
  # Title style (row 1)
  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1, gridExpand = TRUE)

  # Header style (row 2, all columns)
  addStyle(wb, sheet_name, header_style, rows = 2, cols = 1:5, gridExpand = TRUE)

  # Data styles (row 3 onwards)
  n_rows <- nrow(data) + 2  # +2 for title and header rows

  # Text columns (1, 2, 5)
  addStyle(wb, sheet_name, data_style_text, rows = 3:n_rows, cols = text_cols, gridExpand = TRUE)

  # logFC column (3) - 2 decimal places
  addStyle(wb, sheet_name, data_style_logFC, rows = 3:n_rows, cols = 3, gridExpand = TRUE)

  # adj.P.Val column (4) - scientific notation
  addStyle(wb, sheet_name, data_style_pval, rows = 3:n_rows, cols = 4, gridExpand = TRUE)

  # Set column widths
  setColWidths(wb, sheet_name, cols = 1, widths = 16)    # UniProt Accession
  setColWidths(wb, sheet_name, cols = 2, widths = 14)    # Gene Symbol
  setColWidths(wb, sheet_name, cols = 3, widths = 20)    # Avg(log2(Tuni/Ctrl))
  setColWidths(wb, sheet_name, cols = 4, widths = 16)    # Adjusted P value
  setColWidths(wb, sheet_name, cols = 5, widths = 100)   # Protein Annotation

  cat("  Sheet", i, "created with", nrow(data), "proteins\n\n")
}

# =============================================================================
# Save workbook
# =============================================================================

saveWorkbook(wb, output_file, overwrite = TRUE)
cat("Supporting Table S3 saved to:\n", output_file, "\n")
