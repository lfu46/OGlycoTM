# generate_supporting_table_S10_OGalNAc_DE.R
# Supporting Table S10: Abundance changes of O-GalNAcylated proteins

library(tidyverse)
library(openxlsx)

source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'
output_file <- paste0(source_file_path, "supporting_tables/version_Mar24_2026/supporting_table_S10.xlsx")

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
data_style_text <- createStyle(fontName = "Times New Roman", fontSize = 11, halign = "left", valign = "center")
data_style_numeric <- createStyle(fontName = "Times New Roman", fontSize = 11, halign = "center", valign = "center")

wb <- createWorkbook()
cell_types <- c("HEK293T", "HepG2", "Jurkat")
sheet_labels <- c("A", "B", "C")

for (i in seq_along(cell_types)) {
  cell_type <- cell_types[i]
  de <- read_csv(paste0(source_file_path, "differential_analysis/OGalNAc_protein_DE_", cell_type, ".csv"), show_col_types = FALSE)

  result <- de %>%
    dplyr::select(
      `UniProt Accession` = Protein.ID,
      `Gene Symbol` = Gene,
      `Avg(log2(Tuni/Ctrl))` = logFC,
      `Adjusted P value` = adj.P.Val,
      `Protein Annotation` = Protein.Description
    ) %>%
    arrange(`Gene Symbol`)

  sheet_name <- paste0("Sheet", i)
  addWorksheet(wb, sheet_name)
  title_text <- paste0("Table S10", sheet_labels[i],
    ". Abundance changes of O-GalNAcylated proteins in ", cell_type,
    " cells with the inhibition of N-glycosylation")
  writeData(wb, sheet_name, title_text, startRow = 1, startCol = 1, colNames = FALSE)
  mergeCells(wb, sheet_name, cols = 1:5, rows = 1)
  writeData(wb, sheet_name, as.data.frame(t(colnames(result))), startRow = 2, startCol = 1, colNames = FALSE)
  writeData(wb, sheet_name, result, startRow = 3, startCol = 1, colNames = FALSE)

  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1, gridExpand = TRUE)
  addStyle(wb, sheet_name, header_style, rows = 2, cols = 1:5, gridExpand = TRUE)
  n_rows <- nrow(result) + 2
  addStyle(wb, sheet_name, data_style_text, rows = 3:n_rows, cols = c(1,2,5), gridExpand = TRUE)
  addStyle(wb, sheet_name, data_style_numeric, rows = 3:n_rows, cols = c(3,4), gridExpand = TRUE)
  setColWidths(wb, sheet_name, cols = 1, widths = 16)
  setColWidths(wb, sheet_name, cols = 2, widths = 12)
  setColWidths(wb, sheet_name, cols = 3, widths = 20)
  setColWidths(wb, sheet_name, cols = 4, widths = 16)
  setColWidths(wb, sheet_name, cols = 5, widths = 80)

  n_up <- sum(result$`Avg(log2(Tuni/Ctrl))` > 0.5 & result$`Adjusted P value` < 0.05)
  n_down <- sum(result$`Avg(log2(Tuni/Ctrl))` < -0.5 & result$`Adjusted P value` < 0.05)
  cat(cell_type, ":", nrow(result), "proteins,", n_up, "up,", n_down, "down\n")
}

saveWorkbook(wb, output_file, overwrite = TRUE)
cat("Table S10 saved to:", output_file, "\n")
