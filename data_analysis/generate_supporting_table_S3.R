# generate_supporting_table_S3.R
# Generate Supporting Table S3: Identification of O-GlcNAcylation sites
#
# This script creates an Excel file with 3 sheets (HEK293T, HepG2, Jurkat)
# containing O-GlcNAcylation site identification data with proper formatting.

library(tidyverse)
library(openxlsx)

# =============================================================================
# File paths
# =============================================================================

# Input: Site-level O-GlcNAc data
site_files <- c(
  "HEK293T" = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/site/OGlcNAc_site_HEK293T.csv",
  "HepG2" = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/site/OGlcNAc_site_HepG2.csv",
  "Jurkat" = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/site/OGlcNAc_site_Jurkat.csv"
)

# Output file
output_file <- "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/supporting_tables/new_version/supporting_table_S3.xlsx"

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

# =============================================================================
# Process data for each cell type
# =============================================================================

process_cell_type <- function(cell_type, site_file) {

  cat("Processing", cell_type, "...\n")

  # Read site data
  site_data <- read_csv(site_file, show_col_types = FALSE)
  cat("  Total sites:", nrow(site_data), "\n")

  # Select and rename columns to match table format
  result <- site_data %>%
    mutate(
      # Create Site column (e.g., "S402")
      Site = paste0(modified_residue, site_number)
    ) %>%
    select(
      `UniProt Accession` = Protein.ID,
      `Gene Symbol` = Gene,
      Site = Site,
      `Site Position` = site_number,
      `Confidence Level` = Confidence.Level,
      Peptide = Peptide,
      `Obs. m/z` = Observed.M.Z,
      `Precursor Charge` = Charge,
      Hyperscore = Hyperscore,
      `Delta Mass` = Delta.Mass,
      `Total Glycan Composition` = Total.Glycan.Composition,
      `Site Probability` = Site.Probabilities,
      `Protein Annotation` = Protein.Description
    ) %>%
    # Sort by UniProt Accession, then by Site Position
    arrange(`UniProt Accession`, `Site Position`)

  return(result)
}

# =============================================================================
# Create workbook
# =============================================================================

wb <- createWorkbook()

cell_types <- c("HEK293T", "HepG2", "Jurkat")
sheet_labels <- c("A", "B", "C")

# Define column types for formatting
# Text columns: 1 (UniProt), 2 (Gene), 3 (Site), 5 (Confidence Level), 6 (Peptide), 11 (Glycan Comp), 12 (Site Prob), 13 (Annotation)
# Numeric columns: 4 (Site Position), 7 (Obs m/z), 8 (Charge), 9 (Hyperscore), 10 (PPM)
text_cols <- c(1, 2, 3, 5, 6, 11, 12, 13)
numeric_cols <- c(4, 7, 8, 9, 10)

for (i in seq_along(cell_types)) {
  cell_type <- cell_types[i]
  sheet_label <- sheet_labels[i]
  sheet_name <- paste0("Sheet", i)

  # Process data
  data <- process_cell_type(cell_type, site_files[cell_type])

  # Add worksheet
  addWorksheet(wb, sheet_name)

  # Title text
  title_text <- paste0("Table S3", sheet_label, ". Identification of O-GlcNAcylation sites in ", cell_type, " cells")

  # Write title row (row 1)
  writeData(wb, sheet_name, title_text, startRow = 1, startCol = 1, colNames = FALSE)

  # Merge title row across all columns
  mergeCells(wb, sheet_name, cols = 1:13, rows = 1)

  # Write header row (row 2)
  writeData(wb, sheet_name, as.data.frame(t(colnames(data))), startRow = 2, startCol = 1, colNames = FALSE)

  # Write data rows (starting row 3)
  writeData(wb, sheet_name, data, startRow = 3, startCol = 1, colNames = FALSE)

  # Apply styles
  # Title style (row 1)
  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1, gridExpand = TRUE)

  # Header style (row 2, all columns)
  addStyle(wb, sheet_name, header_style, rows = 2, cols = 1:13, gridExpand = TRUE)

  # Data styles (row 3 onwards)
  n_rows <- nrow(data) + 2  # +2 for title and header rows

  # Text columns
  addStyle(wb, sheet_name, data_style_text, rows = 3:n_rows, cols = text_cols, gridExpand = TRUE)

  # Numeric columns
  addStyle(wb, sheet_name, data_style_numeric, rows = 3:n_rows, cols = numeric_cols, gridExpand = TRUE)

  # Set column widths
  setColWidths(wb, sheet_name, cols = 1, widths = 16)    # UniProt Accession
  setColWidths(wb, sheet_name, cols = 2, widths = 12)    # Gene Symbol
  setColWidths(wb, sheet_name, cols = 3, widths = 8)     # Site
  setColWidths(wb, sheet_name, cols = 4, widths = 12)    # Site Position
  setColWidths(wb, sheet_name, cols = 5, widths = 15)    # Confidence Level
  setColWidths(wb, sheet_name, cols = 6, widths = 20)    # Peptide
  setColWidths(wb, sheet_name, cols = 7, widths = 12)    # Obs. m/z
  setColWidths(wb, sheet_name, cols = 8, widths = 15)    # Precursor Charge
  setColWidths(wb, sheet_name, cols = 9, widths = 12)    # Hyperscore
  setColWidths(wb, sheet_name, cols = 10, widths = 10)   # PPM
  setColWidths(wb, sheet_name, cols = 11, widths = 22)   # Total Glycan Composition
  setColWidths(wb, sheet_name, cols = 12, widths = 18)   # Site Probability
  setColWidths(wb, sheet_name, cols = 13, widths = 80)   # Protein Annotation

  cat("  Sheet", i, "created with", nrow(data), "sites\n\n")
}

# =============================================================================
# Save workbook
# =============================================================================

saveWorkbook(wb, output_file, overwrite = TRUE)
cat("Supporting Table S3 saved to:\n", output_file, "\n")
