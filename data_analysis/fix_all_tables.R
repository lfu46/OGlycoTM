# fix_all_tables.R
# Regenerate S1-S7 from source data with proper numeric types
# S8-S11 are already correct

library(tidyverse)
library(openxlsx)

source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'
old_table_dir <- paste0(source_file_path, "supporting_tables/version_Mar24_2026/")
table_dir <- paste0(source_file_path, "supporting_tables/version_Apr5_2026/")
local_dir <- "/Users/longpingfu/Downloads/OGlycoTM/Manuscript/version_Apr5_2026/"

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)

# Copy S8-S11 from old version (not regenerated here)
for (s in 8:11) {
  f <- paste0("supporting_table_S", s, ".xlsx")
  file.copy(paste0(old_table_dir, f), paste0(table_dir, f), overwrite = TRUE)
  file.copy(paste0(old_table_dir, f), paste0(local_dir, f), overwrite = TRUE)
}

cell_types <- c("HEK293T", "HepG2", "Jurkat")
sheet_labels <- c("A", "B", "C")

# =============================================================================
# Shared styles
# =============================================================================
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
data_style_logFC <- createStyle(
  fontName = "Times New Roman", fontSize = 11, halign = "center", valign = "center",
  numFmt = "0.00")
data_style_pval <- createStyle(
  fontName = "Times New Roman", fontSize = 11, halign = "center", valign = "center",
  numFmt = "0.00E+00")
data_style_pct <- createStyle(
  fontName = "Times New Roman", fontSize = 11, halign = "center", valign = "center")
data_style_mz <- createStyle(
  fontName = "Times New Roman", fontSize = 11, halign = "center", valign = "center",
  numFmt = "0.0000")
data_style_hyperscore <- createStyle(
  fontName = "Times New Roman", fontSize = 11, halign = "center", valign = "center",
  numFmt = "0.000")
data_style_deltamass <- createStyle(
  fontName = "Times New Roman", fontSize = 11, halign = "center", valign = "center",
  numFmt = "0.0000")
data_style_percent <- createStyle(
  fontName = "Times New Roman", fontSize = 11, halign = "center", valign = "center",
  numFmt = "0.00%")

# =============================================================================
# Helper: write formatted sheet for protein ID table (8 cols)
# =============================================================================
write_protein_id_sheet <- function(wb, sheet_name, data, title_text, is_wp = FALSE) {
  addWorksheet(wb, sheet_name)
  n_cols <- ncol(data)
  writeData(wb, sheet_name, title_text, startRow = 1, startCol = 1, colNames = FALSE)
  mergeCells(wb, sheet_name, cols = 1:n_cols, rows = 1)
  writeData(wb, sheet_name, as.data.frame(t(colnames(data))), startRow = 2, startCol = 1, colNames = FALSE)
  writeData(wb, sheet_name, data, startRow = 3, startCol = 1, colNames = FALSE)

  n_rows <- nrow(data) + 2
  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1, gridExpand = TRUE)
  addStyle(wb, sheet_name, header_style, rows = 2, cols = 1:n_cols, gridExpand = TRUE)
  addStyle(wb, sheet_name, data_style_text, rows = 3:n_rows, cols = c(1, 2, 8), gridExpand = TRUE)
  addStyle(wb, sheet_name, data_style_numeric, rows = 3:n_rows, cols = 3:6, gridExpand = TRUE)
  addStyle(wb, sheet_name, data_style_percent, rows = 3:n_rows, cols = 7, gridExpand = TRUE)

  setColWidths(wb, sheet_name, cols = 1, widths = 16.3)
  setColWidths(wb, sheet_name, cols = 2, widths = 14.5)
  setColWidths(wb, sheet_name, cols = 3, widths = 11.8)
  setColWidths(wb, sheet_name, cols = 4, widths = 13.6)
  setColWidths(wb, sheet_name, cols = 5, widths = 6.6)
  setColWidths(wb, sheet_name, cols = 6, widths = 8.5)
  setColWidths(wb, sheet_name, cols = 7, widths = 18.3)
  setColWidths(wb, sheet_name, cols = 8, widths = 100)
}

# =============================================================================
# Helper: write formatted sheet for DE table (5 cols)
# =============================================================================
write_de_sheet <- function(wb, sheet_name, data, title_text) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, title_text, startRow = 1, startCol = 1, colNames = FALSE)
  mergeCells(wb, sheet_name, cols = 1:5, rows = 1)
  writeData(wb, sheet_name, as.data.frame(t(colnames(data))), startRow = 2, startCol = 1, colNames = FALSE)
  writeData(wb, sheet_name, data, startRow = 3, startCol = 1, colNames = FALSE)

  n_rows <- nrow(data) + 2
  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1, gridExpand = TRUE)
  addStyle(wb, sheet_name, header_style, rows = 2, cols = 1:5, gridExpand = TRUE)
  addStyle(wb, sheet_name, data_style_text, rows = 3:n_rows, cols = c(1, 2, 5), gridExpand = TRUE)
  addStyle(wb, sheet_name, data_style_logFC, rows = 3:n_rows, cols = 3, gridExpand = TRUE)
  addStyle(wb, sheet_name, data_style_pval, rows = 3:n_rows, cols = 4, gridExpand = TRUE)

  setColWidths(wb, sheet_name, cols = 1, widths = 16)
  setColWidths(wb, sheet_name, cols = 2, widths = 14)
  setColWidths(wb, sheet_name, cols = 3, widths = 20)
  setColWidths(wb, sheet_name, cols = 4, widths = 16)
  setColWidths(wb, sheet_name, cols = 5, widths = 100)
}

# =============================================================================
# Helper: write formatted sheet for site ID table (13 cols)
# =============================================================================
write_site_id_sheet <- function(wb, sheet_name, data, title_text) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, title_text, startRow = 1, startCol = 1, colNames = FALSE)
  mergeCells(wb, sheet_name, cols = 1:13, rows = 1)
  writeData(wb, sheet_name, as.data.frame(t(colnames(data))), startRow = 2, startCol = 1, colNames = FALSE)
  writeData(wb, sheet_name, data, startRow = 3, startCol = 1, colNames = FALSE)

  text_cols <- c(1, 2, 3, 5, 6, 11, 13)
  numeric_cols <- c(4, 7, 8, 9, 10, 12)
  n_rows <- nrow(data) + 2

  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1, gridExpand = TRUE)
  addStyle(wb, sheet_name, header_style, rows = 2, cols = 1:13, gridExpand = TRUE)
  addStyle(wb, sheet_name, data_style_text, rows = 3:n_rows, cols = text_cols, gridExpand = TRUE)
  addStyle(wb, sheet_name, data_style_numeric, rows = 3:n_rows, cols = numeric_cols, gridExpand = TRUE)
  # Override specific numeric columns with fixed decimal formats
  addStyle(wb, sheet_name, data_style_mz, rows = 3:n_rows, cols = 7, gridExpand = TRUE)          # Obs. m/z
  addStyle(wb, sheet_name, data_style_hyperscore, rows = 3:n_rows, cols = 9, gridExpand = TRUE)   # Hyperscore
  addStyle(wb, sheet_name, data_style_deltamass, rows = 3:n_rows, cols = 10, gridExpand = TRUE)   # Delta Mass

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
}

# =============================================================================
# Helper: write formatted sheet for site DE table (6 cols)
# =============================================================================
write_site_de_sheet <- function(wb, sheet_name, data, title_text) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, title_text, startRow = 1, startCol = 1, colNames = FALSE)
  mergeCells(wb, sheet_name, cols = 1:6, rows = 1)
  writeData(wb, sheet_name, as.data.frame(t(colnames(data))), startRow = 2, startCol = 1, colNames = FALSE)
  writeData(wb, sheet_name, data, startRow = 3, startCol = 1, colNames = FALSE)

  n_rows <- nrow(data) + 2
  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1, gridExpand = TRUE)
  addStyle(wb, sheet_name, header_style, rows = 2, cols = 1:6, gridExpand = TRUE)
  addStyle(wb, sheet_name, data_style_text, rows = 3:n_rows, cols = c(1, 2, 3, 6), gridExpand = TRUE)
  addStyle(wb, sheet_name, data_style_logFC, rows = 3:n_rows, cols = 4, gridExpand = TRUE)
  addStyle(wb, sheet_name, data_style_pval, rows = 3:n_rows, cols = 5, gridExpand = TRUE)

  setColWidths(wb, sheet_name, cols = 1, widths = 16)
  setColWidths(wb, sheet_name, cols = 2, widths = 14)
  setColWidths(wb, sheet_name, cols = 3, widths = 10)
  setColWidths(wb, sheet_name, cols = 4, widths = 20)
  setColWidths(wb, sheet_name, cols = 5, widths = 16)
  setColWidths(wb, sheet_name, cols = 6, widths = 100)
}

# =============================================================================
# Helper: save to both locations
# =============================================================================
save_table <- function(wb, table_num) {
  f1 <- paste0(table_dir, "supporting_table_S", table_num, ".xlsx")
  f2 <- paste0(local_dir, "supporting_table_S", table_num, ".xlsx")
  saveWorkbook(wb, f1, overwrite = TRUE)
  file.copy(f1, f2, overwrite = TRUE)
}

# =============================================================================
# S1: Identification of O-GlcNAcylated proteins
# =============================================================================
cat("=== S1: O-GlcNAc protein ID ===\n")

oglcnac_files <- c(
  "HEK293T" = paste0(source_file_path, "filtered/OGlcNAc_HEK293T.csv"),
  "HepG2"   = paste0(source_file_path, "filtered/OGlcNAc_HepG2.csv"),
  "Jurkat"  = paste0(source_file_path, "filtered/OGlcNAc_Jurkat.csv"))

protein_files <- c(
  "HEK293T" = paste0(source_file_path, "../OGlycoTM_HEK293T/OGlyco/EThcD_OPair_TMT_Search/OGlyco/protein.tsv"),
  "HepG2"   = paste0(source_file_path, "../OGlycoTM_HepG2/OGlyco/EThcD_OPair_TMT_Search/OGlyco/protein.tsv"),
  "Jurkat"  = paste0(source_file_path, "../OGlycoTM_Jurkat/OGlyco/EThcD_OPair_TMT_Search/OGlyco/protein.tsv"))

wb <- createWorkbook()
for (i in seq_along(cell_types)) {
  ct <- cell_types[i]
  oglcnac_prots <- read_csv(oglcnac_files[ct], show_col_types = FALSE) |>
    pull(Protein.ID) |> unique()
  prot <- read_tsv(protein_files[ct], show_col_types = FALSE) |>
    filter(`Protein ID` %in% oglcnac_prots) |>
    transmute(
      `UniProt Accession` = `Protein ID`,
      `Gene Symbol` = Gene,
      `Total Peptide` = as.integer(`Total Peptides`),
      `Unique Peptide` = as.integer(`Unique Peptides`),
      Length = as.integer(Length),
      Coverage_pct = Coverage,
      Coverage = as.integer(round(Length * Coverage_pct / 100)),
      `Coverage Percentage` = Coverage_pct / 100,
      `Protein Annotation` = `Protein Description`) |>
    select(-Coverage_pct) |>
    mutate(`Unique Peptide` = ifelse(`Unique Peptide` == 0L, `Total Peptide`, `Unique Peptide`)) |>
    arrange(`UniProt Accession`)
  title <- paste0("Table S1", sheet_labels[i], ". Identification of O-GlcNAcylated proteins in ", ct, " cells")
  write_protein_id_sheet(wb, paste0("Sheet", i), prot, title)
  cat("  ", ct, ":", nrow(prot), "\n")
}
save_table(wb, 1)

# =============================================================================
# S2: Abundance changes of O-GlcNAcylated proteins
# =============================================================================
cat("\n=== S2: O-GlcNAc protein DE ===\n")

wb <- createWorkbook()
for (i in seq_along(cell_types)) {
  ct <- cell_types[i]
  de <- read_csv(paste0(source_file_path, "differential_analysis/OGlcNAc_protein_DE_", ct, ".csv"),
                 show_col_types = FALSE)
  result <- de |>
    transmute(
      `UniProt Accession` = Protein.ID,
      `Gene Symbol` = Gene,
      `Avg(log2(Tuni/Ctrl))` = as.numeric(logFC),
      `Adjusted P value` = as.numeric(adj.P.Val),
      `Protein Annotation` = Protein.Description) |>
    arrange(`UniProt Accession`)
  title <- paste0("Table S2", sheet_labels[i],
    ". Abundance changes of O-GlcNAcylated proteins upon tunicamycin treatment in ", ct, " cells")
  write_de_sheet(wb, paste0("Sheet", i), result, title)
  cat("  ", ct, ":", nrow(result), "\n")
}
save_table(wb, 2)

# =============================================================================
# S3: Identification of total proteins (whole proteome)
# Read from existing table and fix types
# =============================================================================
cat("\n=== S3: WP protein ID ===\n")

wb <- createWorkbook()
for (i in seq_along(cell_types)) {
  ct <- cell_types[i]
  old <- read.xlsx(paste0(old_table_dir, "supporting_table_S3.xlsx"),
                   sheet = i, startRow = 2, colNames = TRUE)
  # Read corrected coverage from FASTA-based computation
  cov_fix <- read_csv(paste0("WP_", ct, "_coverage.csv"), show_col_types = FALSE)
  result <- old |>
    mutate(
      Total.Peptide = as.integer(as.numeric(Total.Peptide)),
      Unique.Peptide = as.integer(as.numeric(Unique.Peptide)),
      Length = as.integer(as.numeric(Length))) |>
    left_join(cov_fix, by = c("UniProt.Accession" = "UniProt_Accession")) |>
    mutate(
      Coverage = as.integer(ifelse(!is.na(Coverage_correct), Coverage_correct, as.numeric(Coverage))),
      Coverage.Percentage = ifelse(Length > 0, Coverage / Length, 0)) |>
    select(-Coverage_correct)
  colnames(result) <- c("UniProt Accession", "Gene Symbol", "Total Peptide",
                        "Unique Peptide", "Length", "Coverage",
                        "Coverage Percentage", "Protein Annotation")

  sheet_name <- paste0("Sheet", i)
  addWorksheet(wb, sheet_name)
  title <- paste0("Table S3", sheet_labels[i], ". Identification of total proteins in ", ct, " cells")
  writeData(wb, sheet_name, title, startRow = 1, startCol = 1, colNames = FALSE)
  mergeCells(wb, sheet_name, cols = 1:8, rows = 1)
  writeData(wb, sheet_name, as.data.frame(t(colnames(result))), startRow = 2, startCol = 1, colNames = FALSE)
  writeData(wb, sheet_name, result, startRow = 3, startCol = 1, colNames = FALSE)

  n_rows <- nrow(result) + 2
  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1, gridExpand = TRUE)
  addStyle(wb, sheet_name, header_style, rows = 2, cols = 1:8, gridExpand = TRUE)
  addStyle(wb, sheet_name, data_style_text, rows = 3:n_rows, cols = c(1, 2, 8), gridExpand = TRUE)
  addStyle(wb, sheet_name, data_style_numeric, rows = 3:n_rows, cols = 3:6, gridExpand = TRUE)
  addStyle(wb, sheet_name, data_style_percent, rows = 3:n_rows, cols = 7, gridExpand = TRUE)

  setColWidths(wb, sheet_name, cols = 1, widths = 16)
  setColWidths(wb, sheet_name, cols = 2, widths = 14)
  setColWidths(wb, sheet_name, cols = 3, widths = 12)
  setColWidths(wb, sheet_name, cols = 4, widths = 14)
  setColWidths(wb, sheet_name, cols = 5, widths = 8)
  setColWidths(wb, sheet_name, cols = 6, widths = 10)
  setColWidths(wb, sheet_name, cols = 7, widths = 18)
  setColWidths(wb, sheet_name, cols = 8, widths = 100)

  cat("  ", ct, ":", nrow(result), "\n")
}
save_table(wb, 3)

# =============================================================================
# S4: Abundance changes of total proteins (whole proteome)
# =============================================================================
cat("\n=== S4: WP protein DE ===\n")

wb <- createWorkbook()
for (i in seq_along(cell_types)) {
  ct <- cell_types[i]
  de <- read_csv(paste0(source_file_path, "differential_analysis/WP_protein_DE_", ct, ".csv"),
                 show_col_types = FALSE)
  result <- de |>
    transmute(
      `UniProt Accession` = UniProt_Accession,
      `Gene Symbol` = Gene.Symbol,
      `Avg(log2(Tuni/Ctrl))` = as.numeric(logFC),
      `Adjusted P value` = as.numeric(adj.P.Val),
      `Protein Annotation` = Annotation) |>
    arrange(`UniProt Accession`)
  title <- paste0("Table S4", sheet_labels[i],
    ". Abundance changes of total proteins upon tunicamycin treatment in ", ct, " cells")
  write_de_sheet(wb, paste0("Sheet", i), result, title)
  cat("  ", ct, ":", nrow(result), "\n")
}
save_table(wb, 4)

# =============================================================================
# S5: Identification of O-GlcNAcylation sites
# =============================================================================
cat("\n=== S5: O-GlcNAc site ID ===\n")

parse_site_probability <- function(prob_str) {
  as.numeric(str_extract(prob_str, "[0-9.]+(?=\\]$)"))
}

wb <- createWorkbook()
for (i in seq_along(cell_types)) {
  ct <- cell_types[i]
  site_data <- read_csv(paste0(source_file_path, "site/OGlcNAc_site_", ct, ".csv"),
                        show_col_types = FALSE) |>
    mutate(Site_Prob = parse_site_probability(Site.Probabilities))

  result <- site_data |>
    filter(Confidence.Level == "Level1" |
           (Confidence.Level == "Level1b" & Site_Prob >= 0.75)) |>
    mutate(Site = paste0(modified_residue, site_number)) |>
    select(
      `UniProt Accession` = Protein.ID,
      `Gene Symbol` = Gene,
      Site,
      `Site Position` = site_number,
      `O-Pair Localization Level` = Confidence.Level,
      Peptide,
      `Obs. m/z` = Observed.M.Z,
      `Precursor Charge` = Charge,
      Hyperscore,
      `Delta Mass` = Delta.Mass,
      `Total Glycan Composition` = Total.Glycan.Composition,
      `Glycosite Position in Peptide` = peptide_site,
      `Protein Annotation` = Protein.Description) |>
    mutate(
      `Site Position` = as.integer(`Site Position`),
      `Obs. m/z` = as.numeric(`Obs. m/z`),
      `Precursor Charge` = as.integer(`Precursor Charge`),
      Hyperscore = as.numeric(Hyperscore),
      `Delta Mass` = as.numeric(`Delta Mass`),
      `Glycosite Position in Peptide` = as.integer(`Glycosite Position in Peptide`)) |>
    arrange(`UniProt Accession`, `Site Position`)

  title <- paste0("Table S5", sheet_labels[i], ". Identification of O-GlcNAcylation sites in ", ct, " cells")
  write_site_id_sheet(wb, paste0("Sheet", i), result, title)
  cat("  ", ct, ":", nrow(result), "\n")
}
save_table(wb, 5)

# =============================================================================
# S6: Abundance changes of O-GlcNAcylation sites
# =============================================================================
cat("\n=== S6: O-GlcNAc site DE ===\n")

wb <- createWorkbook()
for (i in seq_along(cell_types)) {
  ct <- cell_types[i]

  # Get high-confidence site indices
  site_data <- read_csv(paste0(source_file_path, "site/OGlcNAc_site_", ct, ".csv"),
                        show_col_types = FALSE) |>
    mutate(Site_Prob = parse_site_probability(Site.Probabilities))
  high_conf <- site_data |>
    filter(Confidence.Level == "Level1" |
           (Confidence.Level == "Level1b" & Site_Prob >= 0.75)) |>
    pull(site_index) |> unique()

  de <- read_csv(paste0(source_file_path, "differential_analysis/OGlcNAc_site_DE_", ct, ".csv"),
                 show_col_types = FALSE) |>
    filter(site_index %in% high_conf)

  result <- de |>
    mutate(Site = sub("^[^_]+_", "", site_index)) |>
    transmute(
      `UniProt Accession` = Protein.ID,
      `Gene Symbol` = Gene,
      Site,
      `Avg(log2(Tuni/Ctrl))` = as.numeric(logFC),
      `Adjusted P value` = as.numeric(adj.P.Val),
      `Protein Annotation` = Protein.Description) |>
    arrange(`UniProt Accession`, Site)

  title <- paste0("Table S6", sheet_labels[i],
    ". Abundance changes of O-GlcNAcylation sites upon tunicamycin treatment in ", ct, " cells")
  write_site_de_sheet(wb, paste0("Sheet", i), result, title)
  cat("  ", ct, ":", nrow(result), "\n")
}
save_table(wb, 6)

# =============================================================================
# S7: Identification of O-GalNAcylated proteins
# =============================================================================
cat("\n=== S7: O-GalNAc protein ID ===\n")

OGalNAc_compositions <- c(
  "HexNAt(1)GAO_Methoxylamine(1) % 326.1339",
  "HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968")

bonafide_files <- c(
  "HEK293T" = paste0(source_file_path, "filtered/OGlyco_HEK293T_bonafide.csv"),
  "HepG2"   = paste0(source_file_path, "filtered/OGlyco_HepG2_bonafide.csv"),
  "Jurkat"  = paste0(source_file_path, "filtered/OGlyco_Jurkat_bonafide.csv"))

wb <- createWorkbook()
for (i in seq_along(cell_types)) {
  ct <- cell_types[i]
  ogalnac_prots <- read_csv(bonafide_files[ct], show_col_types = FALSE) |>
    filter(Total.Glycan.Composition %in% OGalNAc_compositions) |>
    pull(Protein.ID) |> unique()
  prot <- read_tsv(protein_files[ct], show_col_types = FALSE) |>
    filter(`Protein ID` %in% ogalnac_prots) |>
    transmute(
      `UniProt Accession` = `Protein ID`,
      `Gene Symbol` = Gene,
      `Total Peptide` = as.integer(`Total Peptides`),
      `Unique Peptide` = as.integer(`Unique Peptides`),
      Length = as.integer(Length),
      Coverage_pct = Coverage,
      Coverage = as.integer(round(Length * Coverage_pct / 100)),
      `Coverage Percentage` = Coverage_pct / 100,
      `Protein Annotation` = `Protein Description`) |>
    select(-Coverage_pct) |>
    mutate(`Unique Peptide` = ifelse(`Unique Peptide` == 0L, `Total Peptide`, `Unique Peptide`)) |>
    arrange(`UniProt Accession`)
  title <- paste0("Table S7", sheet_labels[i], ". Identification of O-GalNAcylated proteins in ", ct, " cells")
  write_protein_id_sheet(wb, paste0("Sheet", i), prot, title)
  cat("  ", ct, ":", nrow(prot), "\n")
}
save_table(wb, 7)

# =============================================================================
# S9 & S11: Fix decimal formatting (site ID tables copied from old version)
# =============================================================================
cat("\n=== S9 & S11: Fix decimal formatting ===\n")

for (table_num in c(9, 11)) {
  wb <- loadWorkbook(paste0(local_dir, "supporting_table_S", table_num, ".xlsx"))
  for (i in seq_along(cell_types)) {
    sheet_name <- paste0("Sheet", i)
    d <- read.xlsx(paste0(local_dir, "supporting_table_S", table_num, ".xlsx"),
                   sheet = i, startRow = 2, colNames = TRUE)
    n_rows <- nrow(d) + 2
    # Col 7 = Obs. m/z, Col 9 = Hyperscore, Col 10 = Delta Mass
    addStyle(wb, sheet_name, data_style_mz, rows = 3:n_rows, cols = 7, gridExpand = TRUE)
    addStyle(wb, sheet_name, data_style_hyperscore, rows = 3:n_rows, cols = 9, gridExpand = TRUE)
    addStyle(wb, sheet_name, data_style_deltamass, rows = 3:n_rows, cols = 10, gridExpand = TRUE)
    # Fix S11 title to match manuscript/SI
    if (table_num == 11) {
      new_title <- paste0("Table S11", sheet_labels[i],
        ". Identification of S-GlcNAcylation sites in ", cell_types[i], " cells")
      writeData(wb, sheet_name, new_title, startRow = 1, startCol = 1, colNames = FALSE)
      addStyle(wb, sheet_name, title_style, rows = 1, cols = 1, gridExpand = TRUE)
    }
  }
  saveWorkbook(wb, paste0(table_dir, "supporting_table_S", table_num, ".xlsx"), overwrite = TRUE)
  file.copy(paste0(table_dir, "supporting_table_S", table_num, ".xlsx"),
            paste0(local_dir, "supporting_table_S", table_num, ".xlsx"), overwrite = TRUE)
  cat("  S", table_num, ": formatted\n")
}

# =============================================================================
# Final verification
# =============================================================================
cat("\n=== Verification ===\n")
for (s in 1:11) {
  f <- paste0(local_dir, "supporting_table_S", s, ".xlsx")
  d <- read.xlsx(f, sheet = 1, startRow = 2, colNames = TRUE)
  nc <- ncol(d)
  # Check numeric columns
  num_issues <- 0
  for (j in 1:nc) {
    cls <- class(d[1, j])
    # Columns that should be numeric
    if (nc == 8 && j %in% 3:6) {
      if (cls != "numeric") num_issues <- num_issues + 1
    }
    if (nc == 5 && j %in% 3:4) {
      if (cls != "numeric") num_issues <- num_issues + 1
    }
    if (nc == 13 && j %in% c(4, 7, 8, 9, 10, 12)) {
      if (cls != "numeric") num_issues <- num_issues + 1
    }
    if (nc == 6 && j %in% 4:5) {
      if (cls != "numeric") num_issues <- num_issues + 1
    }
  }
  sizes <- sapply(names(loadWorkbook(f)), function(sn) {
    nrow(read.xlsx(f, sheet = sn, colNames = FALSE)) - 2
  })
  status <- if (num_issues == 0) "OK" else paste0("ISSUES:", num_issues)
  cat(sprintf("  S%-2d: %s  |  %s\n", s,
    paste(paste0(cell_types, "=", sizes), collapse=", "), status))
}
