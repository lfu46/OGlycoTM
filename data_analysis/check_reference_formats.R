# check_reference_formats.R
# Check reference tables for data types and decimal places

library(openxlsx)

ref_dir <- "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/supporting_tables/latest_version/"

for (s in c(1:5, 8)) {
  f <- paste0(ref_dir, "supporting_table_S", s, ".xlsx")
  if (!file.exists(f)) next
  d <- read.xlsx(f, sheet = 1, startRow = 2, colNames = TRUE)
  cat(sprintf("\n=== Reference S%d (%d cols) ===\n", s, ncol(d)))
  for (j in 1:ncol(d)) {
    val <- d[1, j]
    cls <- class(val)
    if (is.numeric(val) && !is.na(val)) {
      # Check multiple rows for decimal consistency
      vals <- d[1:min(5, nrow(d)), j]
      vals <- vals[!is.na(vals)]
      val_strs <- format(vals, scientific = FALSE)
      cat(sprintf("  Col %2d: %-35s class=%-10s samples: %s\n",
        j, colnames(d)[j], cls, paste(val_strs, collapse=", ")))
    } else {
      cat(sprintf("  Col %2d: %-35s class=%-10s sample: %s\n",
        j, colnames(d)[j], cls, as.character(val)))
    }
  }
}

# Also check current tables
cat("\n\n========== CURRENT TABLES ==========\n")
cur_dir <- "/Users/longpingfu/Downloads/OGlycoTM/Manuscript/version_Mar24_2026/"
for (s in 1:11) {
  f <- paste0(cur_dir, "supporting_table_S", s, ".xlsx")
  d <- read.xlsx(f, sheet = 1, startRow = 2, colNames = TRUE)
  cat(sprintf("\n=== Current S%d (%d cols) ===\n", s, ncol(d)))
  for (j in 1:ncol(d)) {
    val <- d[1, j]
    cls <- class(val)
    if (is.numeric(val) && !is.na(val)) {
      vals <- d[1:min(3, nrow(d)), j]
      vals <- vals[!is.na(vals)]
      cat(sprintf("  Col %2d: %-35s class=%-10s samples: %s\n",
        j, colnames(d)[j], cls, paste(vals, collapse=", ")))
    } else {
      cat(sprintf("  Col %2d: %-35s class=%-10s sample: %s\n",
        j, colnames(d)[j], cls, as.character(val)))
    }
  }
}
