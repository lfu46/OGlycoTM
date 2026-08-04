# update_supporting_tables.R
# Renumber existing supporting tables to match revised manuscript numbering
# and generate new O-GalNAc tables (S8, S9, S10)
#
# New numbering:
#   S1:  O-GlcNAc proteins ID         (was S1)
#   S2:  O-GlcNAc protein changes      (was S3)
#   S3:  WP proteins ID                (was S4)
#   S4:  WP protein changes            (was S5)
#   S5:  O-GlcNAc sites ID             (was S6)
#   S6:  O-GlcNAc site changes         (was S7)
#   S7:  O-GalNAc proteins ID          (was S2)
#   S8:  O-GalNAc protein changes      (NEW)
#   S9:  O-GalNAc sites ID             (NEW)
#   S10: O-GalNAc site changes         (NEW)
#   S11: S-GlcNAc (Cys) sites          (was S8)

library(tidyverse)
library(openxlsx)

source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'
table_dir <- paste0(source_file_path, "supporting_tables/version_Mar24_2026/")

# =============================================================================
# Shared styles (same as all existing tables)
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
data_style_logFC <- createStyle(
  fontName = "Times New Roman", fontSize = 11, halign = "center", valign = "center",
  numFmt = "0.00"
)
data_style_pval <- createStyle(
  fontName = "Times New Roman", fontSize = 11, halign = "center", valign = "center",
  numFmt = "0.0E+00"
)

# =============================================================================
# Part 1: Renumber existing tables
# =============================================================================

# Mapping: old_number -> new_number, new_title_base
renumber_map <- list(
  list(old = 1, new = 1, title = "Identification of O-GlcNAcylated proteins in"),
  list(old = 3, new = 2, title = "Abundance changes of O-GlcNAcylated proteins upon tunicamycin treatment in"),
  list(old = 4, new = 3, title = "Identification of total proteins in"),
  list(old = 5, new = 4, title = "Abundance changes of total proteins upon tunicamycin treatment in"),
  list(old = 6, new = 5, title = "Identification of O-GlcNAcylation sites in"),
  list(old = 7, new = 6, title = "Abundance changes of O-GlcNAcylation sites upon tunicamycin treatment in"),
  list(old = 2, new = 7, title = "Identification of O-GalNAcylated proteins in"),
  list(old = 8, new = 11, title = "S-GlcNAcylation sites on cysteine residues identified and removed during data filtering in")
)

cell_types <- c("HEK293T", "HepG2", "Jurkat")
sheet_labels <- c("A", "B", "C")

# Load all old workbooks first (to avoid overwrite conflicts)
old_wbs <- list()
for (entry in renumber_map) {
  old_file <- paste0(table_dir, "supporting_table_S", entry$old, ".xlsx")
  cat("Loading old S", entry$old, "...\n")
  old_wbs[[as.character(entry$old)]] <- loadWorkbook(old_file)
}

# Now read data, create new workbooks with updated titles, and save
for (entry in renumber_map) {
  old_num <- entry$old
  new_num <- entry$new
  title_base <- entry$title

  cat("\nRenumbering S", old_num, " -> S", new_num, "\n")

  old_wb <- old_wbs[[as.character(old_num)]]
  sheet_names <- names(old_wb)

  # Read all data from old workbook
  new_wb <- createWorkbook()

  for (j in seq_along(sheet_names)) {
    sn <- sheet_names[j]
    sheet_label <- sheet_labels[j]

    # Read raw data (no headers, to preserve everything)
    raw <- read.xlsx(old_wb, sheet = sn, colNames = FALSE, skipEmptyRows = FALSE)

    # Determine new title
    new_title <- paste0("Table S", new_num, sheet_label, ". ", title_base, " ", cell_types[j], " cells")

    # Replace title (row 1)
    raw[1, 1] <- new_title

    # Get dimensions
    n_cols <- ncol(raw)
    n_rows <- nrow(raw)

    # Create new sheet
    new_sn <- paste0("Sheet", j)
    addWorksheet(new_wb, new_sn)

    # Write title
    writeData(new_wb, new_sn, new_title, startRow = 1, startCol = 1, colNames = FALSE)
    mergeCells(new_wb, new_sn, cols = 1:n_cols, rows = 1)

    # Write header (row 2)
    writeData(new_wb, new_sn, raw[2, , drop = FALSE], startRow = 2, startCol = 1, colNames = FALSE)

    # Write data (row 3+)
    if (n_rows > 2) {
      writeData(new_wb, new_sn, raw[3:n_rows, , drop = FALSE], startRow = 3, startCol = 1, colNames = FALSE)
    }

    # Apply title and header styles
    addStyle(new_wb, new_sn, title_style, rows = 1, cols = 1, gridExpand = TRUE)
    addStyle(new_wb, new_sn, header_style, rows = 2, cols = 1:n_cols, gridExpand = TRUE)

    # Apply data styles based on table type
    if (n_rows > 2) {
      # Detect table type by number of columns
      if (n_cols == 5) {
        # DE table (UniProt, Gene, logFC, adj.P, Annotation)
        addStyle(new_wb, new_sn, data_style_text, rows = 3:n_rows, cols = c(1, 2, 5), gridExpand = TRUE)
        addStyle(new_wb, new_sn, data_style_logFC, rows = 3:n_rows, cols = 3, gridExpand = TRUE)
        addStyle(new_wb, new_sn, data_style_pval, rows = 3:n_rows, cols = 4, gridExpand = TRUE)

        setColWidths(new_wb, new_sn, cols = 1, widths = 16)
        setColWidths(new_wb, new_sn, cols = 2, widths = 14)
        setColWidths(new_wb, new_sn, cols = 3, widths = 20)
        setColWidths(new_wb, new_sn, cols = 4, widths = 16)
        setColWidths(new_wb, new_sn, cols = 5, widths = 100)

      } else if (n_cols == 6) {
        # Site DE table (UniProt, Gene, Site, logFC, adj.P, Annotation)
        addStyle(new_wb, new_sn, data_style_text, rows = 3:n_rows, cols = c(1, 2, 3, 6), gridExpand = TRUE)
        addStyle(new_wb, new_sn, data_style_logFC, rows = 3:n_rows, cols = 4, gridExpand = TRUE)
        addStyle(new_wb, new_sn, data_style_pval, rows = 3:n_rows, cols = 5, gridExpand = TRUE)

        setColWidths(new_wb, new_sn, cols = 1, widths = 16)
        setColWidths(new_wb, new_sn, cols = 2, widths = 14)
        setColWidths(new_wb, new_sn, cols = 3, widths = 10)
        setColWidths(new_wb, new_sn, cols = 4, widths = 20)
        setColWidths(new_wb, new_sn, cols = 5, widths = 16)
        setColWidths(new_wb, new_sn, cols = 6, widths = 100)

      } else if (n_cols == 8) {
        # Protein ID table (UniProt, Gene, TotalPep, UniquePep, Length, Coverage, Cov%, Annotation)
        addStyle(new_wb, new_sn, data_style_text, rows = 3:n_rows, cols = c(1, 2, 8), gridExpand = TRUE)
        addStyle(new_wb, new_sn, data_style_numeric, rows = 3:n_rows, cols = 3:7, gridExpand = TRUE)

        setColWidths(new_wb, new_sn, cols = 1, widths = 16.3)
        setColWidths(new_wb, new_sn, cols = 2, widths = 14.5)
        setColWidths(new_wb, new_sn, cols = 3, widths = 11.8)
        setColWidths(new_wb, new_sn, cols = 4, widths = 13.6)
        setColWidths(new_wb, new_sn, cols = 5, widths = 6.6)
        setColWidths(new_wb, new_sn, cols = 6, widths = 8.5)
        setColWidths(new_wb, new_sn, cols = 7, widths = 18.3)
        setColWidths(new_wb, new_sn, cols = 8, widths = 100)

      } else if (n_cols == 13) {
        # Site ID table (13 cols)
        text_cols <- c(1, 2, 3, 5, 6, 11, 13)
        numeric_cols <- c(4, 7, 8, 9, 10, 12)
        addStyle(new_wb, new_sn, data_style_text, rows = 3:n_rows, cols = text_cols, gridExpand = TRUE)
        addStyle(new_wb, new_sn, data_style_numeric, rows = 3:n_rows, cols = numeric_cols, gridExpand = TRUE)

        setColWidths(new_wb, new_sn, cols = 1, widths = 16)
        setColWidths(new_wb, new_sn, cols = 2, widths = 12)
        setColWidths(new_wb, new_sn, cols = 3, widths = 8)
        setColWidths(new_wb, new_sn, cols = 4, widths = 12)
        setColWidths(new_wb, new_sn, cols = 5, widths = 15)
        setColWidths(new_wb, new_sn, cols = 6, widths = 20)
        setColWidths(new_wb, new_sn, cols = 7, widths = 12)
        setColWidths(new_wb, new_sn, cols = 8, widths = 15)
        setColWidths(new_wb, new_sn, cols = 9, widths = 12)
        setColWidths(new_wb, new_sn, cols = 10, widths = 10)
        setColWidths(new_wb, new_sn, cols = 11, widths = 22)
        setColWidths(new_wb, new_sn, cols = 12, widths = 18)
        setColWidths(new_wb, new_sn, cols = 13, widths = 80)

      } else {
        # WP or other wide table — apply generic styles
        addStyle(new_wb, new_sn, data_style_text, rows = 3:n_rows, cols = 1:2, gridExpand = TRUE)
        addStyle(new_wb, new_sn, data_style_numeric, rows = 3:n_rows, cols = 3:n_cols, gridExpand = TRUE)
      }
    }

    cat("  Sheet", j, "(", cell_types[j], "):", n_rows - 2, "data rows\n")
  }

  new_file <- paste0(table_dir, "supporting_table_S", new_num, ".xlsx")
  saveWorkbook(new_wb, new_file, overwrite = TRUE)
  cat("  Saved:", new_file, "\n")
}

# =============================================================================
# Part 2: Generate new S8 — O-GalNAc protein abundance changes
# =============================================================================

cat("\n=== Generating new S8: O-GalNAc protein abundance changes ===\n")

wb_s8 <- createWorkbook()

for (i in seq_along(cell_types)) {
  cell_type <- cell_types[i]
  de_file <- paste0(source_file_path, "differential_analysis/OGalNAc_protein_DE_", cell_type, ".csv")
  de <- read_csv(de_file, show_col_types = FALSE)

  result <- de %>%
    select(
      `UniProt Accession` = Protein.ID,
      `Gene Symbol` = Gene,
      `Avg(log2(Tuni/Ctrl))` = logFC,
      `Adjusted P value` = adj.P.Val,
      `Protein Annotation` = Protein.Description
    ) %>%
    arrange(`UniProt Accession`)

  sheet_name <- paste0("Sheet", i)
  addWorksheet(wb_s8, sheet_name)
  title_text <- paste0("Table S8", sheet_labels[i],
    ". Abundance changes of O-GalNAcylated proteins upon tunicamycin treatment in ",
    cell_type, " cells")

  writeData(wb_s8, sheet_name, title_text, startRow = 1, startCol = 1, colNames = FALSE)
  mergeCells(wb_s8, sheet_name, cols = 1:5, rows = 1)
  writeData(wb_s8, sheet_name, as.data.frame(t(colnames(result))), startRow = 2, startCol = 1, colNames = FALSE)
  writeData(wb_s8, sheet_name, result, startRow = 3, startCol = 1, colNames = FALSE)

  n_rows <- nrow(result) + 2
  addStyle(wb_s8, sheet_name, title_style, rows = 1, cols = 1, gridExpand = TRUE)
  addStyle(wb_s8, sheet_name, header_style, rows = 2, cols = 1:5, gridExpand = TRUE)
  addStyle(wb_s8, sheet_name, data_style_text, rows = 3:n_rows, cols = c(1, 2, 5), gridExpand = TRUE)
  addStyle(wb_s8, sheet_name, data_style_logFC, rows = 3:n_rows, cols = 3, gridExpand = TRUE)
  addStyle(wb_s8, sheet_name, data_style_pval, rows = 3:n_rows, cols = 4, gridExpand = TRUE)

  setColWidths(wb_s8, sheet_name, cols = 1, widths = 16)
  setColWidths(wb_s8, sheet_name, cols = 2, widths = 14)
  setColWidths(wb_s8, sheet_name, cols = 3, widths = 20)
  setColWidths(wb_s8, sheet_name, cols = 4, widths = 16)
  setColWidths(wb_s8, sheet_name, cols = 5, widths = 100)

  cat("  S8", sheet_labels[i], "(", cell_type, "):", nrow(result), "proteins\n")
}

saveWorkbook(wb_s8, paste0(table_dir, "supporting_table_S8.xlsx"), overwrite = TRUE)
cat("  Saved S8\n")

# =============================================================================
# Part 3: Generate new S9 — O-GalNAc site identification
# =============================================================================

cat("\n=== Generating new S9: O-GalNAc site identification ===\n")

parse_site_probability <- function(prob_str) {
  prob <- str_extract(prob_str, "[0-9.]+(?=\\]$)")
  as.numeric(prob)
}

wb_s9 <- createWorkbook()

for (i in seq_along(cell_types)) {
  cell_type <- cell_types[i]
  site_file <- paste0(source_file_path, "site/OGalNAc_site_", cell_type, ".csv")
  site_data <- read_csv(site_file, show_col_types = FALSE)
  cat("  ", cell_type, "raw PSMs:", nrow(site_data), "\n")

  site_data <- site_data %>%
    mutate(Site_Prob = parse_site_probability(Site.Probabilities))

  site_data_filtered <- site_data %>%
    filter(
      Confidence.Level == "Level1" |
      (Confidence.Level == "Level1b" & Site_Prob >= 0.75)
    )

  cat("    After filtering:", nrow(site_data_filtered), "PSMs\n")

  result <- site_data_filtered %>%
    mutate(Site = paste0(modified_residue, site_number)) %>%
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
      `Site Probability` = Site_Prob,
      `Protein Annotation` = Protein.Description
    ) %>%
    arrange(`UniProt Accession`, `Site Position`)

  sheet_name <- paste0("Sheet", i)
  addWorksheet(wb_s9, sheet_name)
  title_text <- paste0("Table S9", sheet_labels[i],
    ". Identification of O-GalNAcylation sites in ", cell_type, " cells")

  writeData(wb_s9, sheet_name, title_text, startRow = 1, startCol = 1, colNames = FALSE)
  mergeCells(wb_s9, sheet_name, cols = 1:13, rows = 1)
  writeData(wb_s9, sheet_name, as.data.frame(t(colnames(result))), startRow = 2, startCol = 1, colNames = FALSE)
  writeData(wb_s9, sheet_name, result, startRow = 3, startCol = 1, colNames = FALSE)

  n_rows <- nrow(result) + 2
  text_cols <- c(1, 2, 3, 5, 6, 11, 13)
  numeric_cols <- c(4, 7, 8, 9, 10, 12)

  addStyle(wb_s9, sheet_name, title_style, rows = 1, cols = 1, gridExpand = TRUE)
  addStyle(wb_s9, sheet_name, header_style, rows = 2, cols = 1:13, gridExpand = TRUE)
  addStyle(wb_s9, sheet_name, data_style_text, rows = 3:n_rows, cols = text_cols, gridExpand = TRUE)
  addStyle(wb_s9, sheet_name, data_style_numeric, rows = 3:n_rows, cols = numeric_cols, gridExpand = TRUE)

  setColWidths(wb_s9, sheet_name, cols = 1, widths = 16)
  setColWidths(wb_s9, sheet_name, cols = 2, widths = 12)
  setColWidths(wb_s9, sheet_name, cols = 3, widths = 8)
  setColWidths(wb_s9, sheet_name, cols = 4, widths = 12)
  setColWidths(wb_s9, sheet_name, cols = 5, widths = 15)
  setColWidths(wb_s9, sheet_name, cols = 6, widths = 20)
  setColWidths(wb_s9, sheet_name, cols = 7, widths = 12)
  setColWidths(wb_s9, sheet_name, cols = 8, widths = 15)
  setColWidths(wb_s9, sheet_name, cols = 9, widths = 12)
  setColWidths(wb_s9, sheet_name, cols = 10, widths = 10)
  setColWidths(wb_s9, sheet_name, cols = 11, widths = 22)
  setColWidths(wb_s9, sheet_name, cols = 12, widths = 18)
  setColWidths(wb_s9, sheet_name, cols = 13, widths = 80)

  cat("    S9", sheet_labels[i], ":", nrow(result), "sites\n")
}

saveWorkbook(wb_s9, paste0(table_dir, "supporting_table_S9.xlsx"), overwrite = TRUE)
cat("  Saved S9\n")

# =============================================================================
# Part 4: Generate new S10 — O-GalNAc site abundance changes
# =============================================================================

cat("\n=== Generating new S10: O-GalNAc site abundance changes ===\n")

wb_s10 <- createWorkbook()

for (i in seq_along(cell_types)) {
  cell_type <- cell_types[i]

  # Read site DE data
  de_file <- paste0(source_file_path, "differential_analysis/OGalNAc_site_DE_", cell_type, ".csv")
  de_data <- read_csv(de_file, show_col_types = FALSE)

  # Read site PSM data for probability filtering
  site_file <- paste0(source_file_path, "site/OGalNAc_site_", cell_type, ".csv")
  site_data <- read_csv(site_file, show_col_types = FALSE) %>%
    mutate(Site_Prob = parse_site_probability(Site.Probabilities))

  # Get high-confidence site indices
  high_conf_sites <- site_data %>%
    filter(
      Confidence.Level == "Level1" |
      (Confidence.Level == "Level1b" & Site_Prob >= 0.75)
    ) %>%
    pull(site_index) %>%
    unique()

  # Filter DE data to high-confidence sites
  de_filtered <- de_data %>%
    filter(site_index %in% high_conf_sites)

  cat("  ", cell_type, ": ", nrow(de_data), " -> ", nrow(de_filtered), " sites after filtering\n")

  result <- de_filtered %>%
    mutate(Site = sub("^[^_]+_", "", site_index)) %>%
    select(
      `UniProt Accession` = Protein.ID,
      `Gene Symbol` = Gene,
      Site = Site,
      `Avg(log2(Tuni/Ctrl))` = logFC,
      `Adjusted P value` = adj.P.Val,
      `Protein Annotation` = Protein.Description
    ) %>%
    arrange(`UniProt Accession`, Site)

  sheet_name <- paste0("Sheet", i)
  addWorksheet(wb_s10, sheet_name)
  title_text <- paste0("Table S10", sheet_labels[i],
    ". Abundance changes of O-GalNAcylation sites upon tunicamycin treatment in ",
    cell_type, " cells")

  writeData(wb_s10, sheet_name, title_text, startRow = 1, startCol = 1, colNames = FALSE)
  mergeCells(wb_s10, sheet_name, cols = 1:6, rows = 1)
  writeData(wb_s10, sheet_name, as.data.frame(t(colnames(result))), startRow = 2, startCol = 1, colNames = FALSE)
  writeData(wb_s10, sheet_name, result, startRow = 3, startCol = 1, colNames = FALSE)

  n_rows <- nrow(result) + 2
  addStyle(wb_s10, sheet_name, title_style, rows = 1, cols = 1, gridExpand = TRUE)
  addStyle(wb_s10, sheet_name, header_style, rows = 2, cols = 1:6, gridExpand = TRUE)
  addStyle(wb_s10, sheet_name, data_style_text, rows = 3:n_rows, cols = c(1, 2, 3, 6), gridExpand = TRUE)
  addStyle(wb_s10, sheet_name, data_style_logFC, rows = 3:n_rows, cols = 4, gridExpand = TRUE)
  addStyle(wb_s10, sheet_name, data_style_pval, rows = 3:n_rows, cols = 5, gridExpand = TRUE)

  setColWidths(wb_s10, sheet_name, cols = 1, widths = 16)
  setColWidths(wb_s10, sheet_name, cols = 2, widths = 14)
  setColWidths(wb_s10, sheet_name, cols = 3, widths = 10)
  setColWidths(wb_s10, sheet_name, cols = 4, widths = 20)
  setColWidths(wb_s10, sheet_name, cols = 5, widths = 16)
  setColWidths(wb_s10, sheet_name, cols = 6, widths = 100)

  cat("    S10", sheet_labels[i], ":", nrow(result), "sites\n")
}

saveWorkbook(wb_s10, paste0(table_dir, "supporting_table_S10.xlsx"), overwrite = TRUE)
cat("  Saved S10\n")

# =============================================================================
# Final summary
# =============================================================================

cat("\n=== All supporting tables updated ===\n")
cat("Output directory:", table_dir, "\n")
for (s in 1:11) {
  f <- paste0(table_dir, "supporting_table_S", s, ".xlsx")
  if (file.exists(f)) {
    wb_check <- loadWorkbook(f)
    sizes <- sapply(names(wb_check), function(sn) {
      d <- read.xlsx(wb_check, sheet = sn, colNames = FALSE)
      nrow(d) - 2  # subtract title + header
    })
    cat(sprintf("  S%-2d: %s\n", s, paste(paste0(cell_types, "=", sizes), collapse = ", ")))
  }
}
