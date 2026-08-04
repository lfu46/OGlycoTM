# regenerate_S11.R
# Regenerate Table S11 from source data with proper numeric types
# Matches S5 format (column headers, data types)

library(tidyverse)
library(openxlsx)

source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'
table_dir <- paste0(source_file_path, "supporting_tables/version_Mar24_2026/")
local_dir <- "/Users/longpingfu/Downloads/OGlycoTM/Manuscript/version_Mar24_2026/"

# Styles
title_style <- createStyle(
  fontName = "Times New Roman", fontSize = 16, fontColour = "#0070C0",
  textDecoration = "bold", halign = "left", valign = "center")
header_style <- createStyle(
  fontName = "Times New Roman", fontSize = 12, fontColour = "#000000",
  textDecoration = "bold", halign = "center", valign = "center",
  border = c("top", "bottom"), borderStyle = c("thin", "double"),
  borderColour = c("#000000", "#000000"))
data_style_text <- createStyle(
  fontName = "Times New Roman", fontSize = 11, halign = "left", valign = "center")
data_style_numeric <- createStyle(
  fontName = "Times New Roman", fontSize = 11, halign = "center", valign = "center")

raw_files <- c(
  "HEK293T" = paste0(source_file_path, "raw/OGlyco_HEK293T_raw.csv"),
  "HepG2"   = paste0(source_file_path, "raw/OGlyco_HepG2_raw.csv"),
  "Jurkat"  = paste0(source_file_path, "raw/OGlyco_Jurkat_raw.csv")
)

cell_types <- c("HEK293T", "HepG2", "Jurkat")
sheet_labels <- c("A", "B", "C")
text_cols <- c(1, 2, 3, 5, 6, 11, 13)
numeric_cols <- c(4, 7, 8, 9, 10, 12)

wb <- createWorkbook()

for (i in seq_along(cell_types)) {
  cell_type <- cell_types[i]
  raw_data <- read_csv(raw_files[cell_type], show_col_types = FALSE)

  # Level1/1b only, Cys sites with HexNAc modification
  cys_pattern_528 <- "C\\(528\\.2859\\)"
  cys_pattern_299 <- "C\\(299\\.1230\\)"

  cys_sglcnac <- raw_data |>
    filter(Confidence.Level %in% c("Level1", "Level1b")) |>
    filter(str_detect(Assigned.Modifications, cys_pattern_528) |
           str_detect(Assigned.Modifications, cys_pattern_299))

  # Extract Cys modification position
  cys_mod_pattern <- "\\d+C\\((528\\.2859|299\\.1230)\\)"

  result <- cys_sglcnac |>
    mutate(
      Cys_Pep_Pos = as.integer(str_extract(
        str_extract(Assigned.Modifications, cys_mod_pattern), "^\\d+")),
      Site_Position = as.integer(Protein.Start + Cys_Pep_Pos - 1),
      Site = paste0("C", Site_Position)
    ) |>
    select(
      `UniProt Accession` = Protein.ID,
      `Gene Symbol` = Gene,
      Site,
      `Site Position` = Site_Position,
      `O-Pair Localization Level` = Confidence.Level,
      Peptide,
      `Obs. m/z` = Observed.M.Z,
      `Precursor Charge` = Charge,
      Hyperscore,
      `Delta Mass` = Delta.Mass,
      `Total Glycan Composition` = Total.Glycan.Composition,
      `Glycosite Position in Peptide` = Cys_Pep_Pos,
      `Protein Annotation` = Protein.Description
    ) |>
    # Best scoring PSM per unique site
    group_by(`UniProt Accession`, Peptide, `Site Position`) |>
    slice_max(Hyperscore, n = 1, with_ties = FALSE) |>
    ungroup() |>
    # Ensure proper numeric types
    mutate(
      `Site Position` = as.integer(`Site Position`),
      `Obs. m/z` = as.numeric(`Obs. m/z`),
      `Precursor Charge` = as.integer(`Precursor Charge`),
      Hyperscore = as.numeric(Hyperscore),
      `Delta Mass` = as.numeric(`Delta Mass`),
      `Glycosite Position in Peptide` = as.integer(`Glycosite Position in Peptide`)
    ) |>
    arrange(`Gene Symbol`, `Site Position`)

  sheet_name <- paste0("Sheet", i)
  addWorksheet(wb, sheet_name)
  title_text <- paste0("Table S11", sheet_labels[i],
    ". S-GlcNAcylation sites on cysteine residues identified and removed during data filtering in ",
    cell_type, " cells")
  writeData(wb, sheet_name, title_text, startRow = 1, startCol = 1, colNames = FALSE)
  mergeCells(wb, sheet_name, cols = 1:13, rows = 1)
  writeData(wb, sheet_name, as.data.frame(t(colnames(result))), startRow = 2, startCol = 1, colNames = FALSE)
  writeData(wb, sheet_name, result, startRow = 3, startCol = 1, colNames = FALSE)

  n_rows <- nrow(result) + 2
  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1, gridExpand = TRUE)
  addStyle(wb, sheet_name, header_style, rows = 2, cols = 1:13, gridExpand = TRUE)
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

  cat(cell_type, ":", nrow(result), "sites\n")
}

saveWorkbook(wb, paste0(table_dir, "supporting_table_S11.xlsx"), overwrite = TRUE)
file.copy(paste0(table_dir, "supporting_table_S11.xlsx"),
          paste0(local_dir, "supporting_table_S11.xlsx"), overwrite = TRUE)
cat("S11 saved to both locations\n")
