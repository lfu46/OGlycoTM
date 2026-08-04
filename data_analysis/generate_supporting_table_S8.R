# generate_supporting_table_S8.R
# Generate Supporting Table S8: S-GlcNAcylation sites on cysteine residues
# removed during data filtering
#
# These are Level1/1b glycopeptides with glycan localized on cysteine,
# which were excluded from the bonafide O-GlcNAc dataset.

library(tidyverse)
library(openxlsx)

# =============================================================================
# File paths
# =============================================================================

source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'

raw_files <- c(
  "HEK293T" = paste0(source_file_path, "raw/OGlyco_HEK293T_raw.csv"),
  "HepG2" = paste0(source_file_path, "raw/OGlyco_HepG2_raw.csv"),
  "Jurkat" = paste0(source_file_path, "raw/OGlyco_Jurkat_raw.csv")
)

output_file <- paste0(source_file_path, "supporting_tables/version_Mar24_2026/supporting_table_S8.xlsx")

# =============================================================================
# Define styles (same as S1-S6)
# =============================================================================

title_style <- createStyle(
  fontName = "Times New Roman",
  fontSize = 16,
  fontColour = "#0070C0",
  textDecoration = "bold",
  halign = "left",
  valign = "center"
)

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

data_style_text <- createStyle(
  fontName = "Times New Roman",
  fontSize = 11,
  halign = "left",
  valign = "center"
)

data_style_numeric <- createStyle(
  fontName = "Times New Roman",
  fontSize = 11,
  halign = "center",
  valign = "center"
)

# =============================================================================
# Process data for each cell type
# =============================================================================

process_cell_type <- function(cell_type, raw_file) {

  cat("Processing", cell_type, "...\n")

  raw_data <- read_csv(raw_file, show_col_types = FALSE)
  cat("  Total raw PSMs:", nrow(raw_data), "\n")

  # Level1/1b only
  localized <- raw_data |>
    filter(Confidence.Level %in% c("Level1", "Level1b"))
  cat("  Level1/1b PSMs:", nrow(localized), "\n")

  # S-GlcNAc on Cys: C(528.2859) or C(299.1230)
  cys_sglcnac <- localized |>
    filter(
      str_detect(Assigned.Modifications, "C\\(528\\.2859\\)") |
      str_detect(Assigned.Modifications, "C\\(299\\.1230\\)")
    )
  cat("  S-GlcNAc Cys PSMs:", nrow(cys_sglcnac), "\n")

  # Extract Cys site positions
  result <- cys_sglcnac |>
    mutate(
      # Extract position of Cys modification
      Cys_Pep_Pos = as.integer(str_extract(
        str_extract(Assigned.Modifications, "\\d+C\\((528\\.2859|299\\.1230)\\)"),
        "^\\d+"
      )),
      # Calculate protein position
      Site_Position = Protein.Start + Cys_Pep_Pos - 1,
      # Create Site label
      Site = paste0("C", Site_Position)
    ) |>
    select(
      `UniProt Accession` = Protein.ID,
      `Gene Symbol` = Gene,
      Site = Site,
      `Site Position` = Site_Position,
      `Confidence Level` = Confidence.Level,
      Peptide = Peptide,
      `Obs. m/z` = Observed.M.Z,
      `Precursor Charge` = Charge,
      Hyperscore = Hyperscore,
      `Delta Mass` = Delta.Mass,
      `Total Glycan Composition` = Total.Glycan.Composition,
      `Glycosite Position in Peptide` = Cys_Pep_Pos,
      `Protein Annotation` = Protein.Description
    ) |>
    # Keep best scoring PSM per unique site
    group_by(`UniProt Accession`, Peptide, `Site Position`) |>
    slice_max(Hyperscore, n = 1, with_ties = FALSE) |>
    ungroup() |>
    arrange(`Gene Symbol`, `Site Position`)

  cat("  Unique Cys sites:", nrow(result), "\n\n")
  return(result)
}

# =============================================================================
# Create workbook
# =============================================================================

wb <- createWorkbook()

cell_types <- c("HEK293T", "HepG2", "Jurkat")
sheet_labels <- c("A", "B", "C")

text_cols <- c(1, 2, 3, 5, 6, 11, 13)
numeric_cols <- c(4, 7, 8, 9, 10, 12)

for (i in seq_along(cell_types)) {
  cell_type <- cell_types[i]
  sheet_label <- sheet_labels[i]
  sheet_name <- paste0("Sheet", i)

  data <- process_cell_type(cell_type, raw_files[cell_type])

  addWorksheet(wb, sheet_name)

  title_text <- paste0(
    "Table S8", sheet_label,
    ". S-GlcNAcylation sites on cysteine residues identified and removed during data filtering in ",
    cell_type, " cells"
  )

  writeData(wb, sheet_name, title_text, startRow = 1, startCol = 1, colNames = FALSE)
  mergeCells(wb, sheet_name, cols = 1:13, rows = 1)
  writeData(wb, sheet_name, as.data.frame(t(colnames(data))), startRow = 2, startCol = 1, colNames = FALSE)
  writeData(wb, sheet_name, data, startRow = 3, startCol = 1, colNames = FALSE)

  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1, gridExpand = TRUE)
  addStyle(wb, sheet_name, header_style, rows = 2, cols = 1:13, gridExpand = TRUE)

  n_rows <- nrow(data) + 2
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

  cat("  Sheet", i, "created with", nrow(data), "sites\n\n")
}

saveWorkbook(wb, output_file, overwrite = TRUE)
cat("Supporting Table S8 saved to:\n", output_file, "\n")
