# generate_supporting_table_S9_OGalNAc.R
# Generate Supporting Table S9: O-GalNAc site identification and quantification
# Format matches Table S6 (O-GlcNAc sites)

library(tidyverse)
library(openxlsx)

# =============================================================================
# File paths
# =============================================================================
source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'

site_files <- c(
  "HEK293T" = paste0(source_file_path, "site/OGalNAc_site_HEK293T.csv"),
  "HepG2" = paste0(source_file_path, "site/OGalNAc_site_HepG2.csv"),
  "Jurkat" = paste0(source_file_path, "site/OGalNAc_site_Jurkat.csv")
)

output_file <- paste0(source_file_path, "supporting_tables/version_Mar24_2026/supporting_table_S9.xlsx")

# =============================================================================
# Styles (same as S6)
# =============================================================================
title_style <- createStyle(
  fontName = "Times New Roman", fontSize = 16, fontColour = "#0070C0",
  textDecoration = "bold", halign = "left", valign = "center"
)
header_style <- createStyle(
  fontName = "Times New Roman", fontSize = 12, fontColour = "#000000",
  textDecoration = "bold", halign = "center", valign = "center",
  border = c("top", "bottom"), borderStyle = c("thin", "double"),
  borderColour = c("#000000", "#000000")
)
data_style_text <- createStyle(
  fontName = "Times New Roman", fontSize = 11, halign = "left", valign = "center"
)
data_style_numeric <- createStyle(
  fontName = "Times New Roman", fontSize = 11, halign = "center", valign = "center"
)

# =============================================================================
# Helper
# =============================================================================
parse_site_probability <- function(prob_str) {
  prob <- str_extract(prob_str, "[0-9.]+(?=\\]$)")
  as.numeric(prob)
}

# =============================================================================
# Process
# =============================================================================
wb <- createWorkbook()
cell_types <- c("HEK293T", "HepG2", "Jurkat")
sheet_labels <- c("A", "B", "C")
text_cols <- c(1, 2, 3, 5, 6, 11, 13)
numeric_cols <- c(4, 7, 8, 9, 10, 12)

for (i in seq_along(cell_types)) {
  cell_type <- cell_types[i]
  sheet_label <- sheet_labels[i]

  site_data <- read_csv(site_files[cell_type], show_col_types = FALSE)
  cat("Processing", cell_type, ":", nrow(site_data), "PSMs\n")

  # Parse probability and filter (Level1 all, Level1b >= 0.75)
  site_data <- site_data %>%
    mutate(Site_Prob = parse_site_probability(Site.Probabilities))

  site_data_filtered <- site_data %>%
    filter(
      Confidence.Level == "Level1" |
      (Confidence.Level == "Level1b" & Site_Prob >= 0.75)
    )

  cat("  After filtering:", nrow(site_data_filtered), "PSMs\n")

  result <- site_data_filtered %>%
    mutate(Site = paste0(modified_residue, site_number)) %>%
    dplyr::select(
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
      `Site Probability` = Site_Prob,
      `Protein Annotation` = Protein.Description
    ) %>%
    arrange(`UniProt Accession`, `Site Position`)

  # Write sheet
  sheet_name <- paste0("Sheet", i)
  addWorksheet(wb, sheet_name)
  title_text <- paste0("Table S9", sheet_label, ". Identification of O-GalNAcylation sites in ", cell_type, " cells")
  writeData(wb, sheet_name, title_text, startRow = 1, startCol = 1, colNames = FALSE)
  mergeCells(wb, sheet_name, cols = 1:13, rows = 1)
  writeData(wb, sheet_name, as.data.frame(t(colnames(result))), startRow = 2, startCol = 1, colNames = FALSE)
  writeData(wb, sheet_name, result, startRow = 3, startCol = 1, colNames = FALSE)

  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1, gridExpand = TRUE)
  addStyle(wb, sheet_name, header_style, rows = 2, cols = 1:13, gridExpand = TRUE)
  n_rows <- nrow(result) + 2
  addStyle(wb, sheet_name, data_style_text, rows = 3:n_rows, cols = text_cols, gridExpand = TRUE)
  addStyle(wb, sheet_name, data_style_numeric, rows = 3:n_rows, cols = numeric_cols, gridExpand = TRUE)

  setColWidths(wb, sheet_name, cols = 1, widths = 16)
  setColWidths(wb, sheet_name, cols = 2, widths = 12)
  setColWidths(wb, sheet_name, cols = 3, widths = 8)
  setColWidths(wb, sheet_name, cols = 4, widths = 12)
  setColWidths(wb, sheet_name, cols = 5, widths = 15)
  setColWidths(wb, sheet_name, cols = 6, widths = 20)
  setColWidths(wb, sheet_name, cols = 7, widths = 12)
  setColWidths(wb, sheet_name, cols = 8, widths = 15)
  setColWidths(wb, sheet_name, cols = 9, widths = 12)
  setColWidths(wb, sheet_name, cols = 10, widths = 10)
  setColWidths(wb, sheet_name, cols = 11, widths = 22)
  setColWidths(wb, sheet_name, cols = 12, widths = 18)
  setColWidths(wb, sheet_name, cols = 13, widths = 80)

  cat("  Sheet", i, ":", nrow(result), "sites\n\n")
}

saveWorkbook(wb, output_file, overwrite = TRUE)
cat("Table S9 saved to:", output_file, "\n")
