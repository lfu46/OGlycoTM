# Find alternative examples
# Criteria:
# - IDR site: logFC slightly below 0 (-0.3 to 0)
# - Structured site: LARGEST positive logFC (any value)

library(tidyverse)
source('data_source.R')

features_analysis <- read_csv(
  paste0(source_file_path, 'site_features/OGlcNAc_site_features.csv'),
  show_col_types = FALSE
)

cat("=== SEARCH: IDR site (-0.3 to 0) + Structured site (max positive) ===\n\n")

for (cell_type in c("HEK293T", "HepG2", "Jurkat")) {

  cat(rep("=", 60), "\n", sep = "")
  cat("Cell Type:", cell_type, "\n")
  cat(rep("=", 60), "\n\n")

  cell_features <- features_analysis |>
    filter(cell == cell_type, !is.na(is_IDR), !is.na(logFC))

  # IDR sites with logFC slightly below 0
  idr_ideal <- cell_features |>
    filter(is_IDR == 1, logFC >= -0.3, logFC <= 0) |>
    select(Protein.ID, Gene, site_number, logFC, pLDDT)

  # Structured sites with positive logFC (any positive value)
  structured_positive <- cell_features |>
    filter(is_IDR == 0, logFC > 0) |>
    select(Protein.ID, Gene, site_number, logFC, pLDDT) |>
    arrange(desc(logFC))

  # Find common proteins (excluding TAB2)
  common_proteins <- setdiff(
    intersect(idr_ideal$Protein.ID, structured_positive$Protein.ID),
    "Q9NYJ8"
  )

  if (length(common_proteins) > 0) {
    cat("Found", length(common_proteins), "proteins (excluding TAB2):\n\n")

    for (pid in common_proteins) {
      gene <- cell_features |> filter(Protein.ID == pid) |> pull(Gene) |> unique()

      # Get the structured site info
      struct_info <- structured_positive |> filter(Protein.ID == pid)
      idr_info <- idr_ideal |> filter(Protein.ID == pid)

      cat(">>> ", pid, " (", gene, ") <<<\n", sep = "")
      cat("Structured site: S", struct_info$site_number[1],
          " logFC = ", round(struct_info$logFC[1], 3),
          " pLDDT = ", round(struct_info$pLDDT[1], 1), "\n", sep = "")
      cat("IDR site: S", idr_info$site_number[1],
          " logFC = ", round(idr_info$logFC[1], 3),
          " pLDDT = ", round(idr_info$pLDDT[1], 1), "\n\n", sep = "")
    }
  } else {
    cat("No matching proteins (excluding TAB2)\n\n")
  }
}

# ============================================
# EXPANDED SEARCH: IDR logFC -0.5 to 0.1
# ============================================

cat("\n")
cat(rep("=", 60), "\n", sep = "")
cat("=== EXPANDED: IDR site (-0.5 to 0.1) + Structured site (positive) ===\n")
cat(rep("=", 60), "\n\n")

for (cell_type in c("HEK293T", "HepG2", "Jurkat")) {

  cat("--- ", cell_type, " ---\n", sep = "")

  cell_features <- features_analysis |>
    filter(cell == cell_type, !is.na(is_IDR), !is.na(logFC))

  # Expanded IDR criteria
  idr_expanded <- cell_features |>
    filter(is_IDR == 1, logFC >= -0.5, logFC <= 0.1) |>
    select(Protein.ID, Gene, site_number, logFC, pLDDT)

  # Structured sites with positive logFC
  structured_positive <- cell_features |>
    filter(is_IDR == 0, logFC > 0) |>
    select(Protein.ID, Gene, site_number, logFC, pLDDT) |>
    arrange(desc(logFC))

  # Find common proteins (excluding TAB2)
  common_proteins <- setdiff(
    intersect(idr_expanded$Protein.ID, structured_positive$Protein.ID),
    "Q9NYJ8"
  )

  if (length(common_proteins) > 0) {
    cat("Found", length(common_proteins), "proteins:\n\n")

    # Sort by structured logFC (largest first)
    results <- tibble()
    for (pid in common_proteins) {
      struct_info <- structured_positive |> filter(Protein.ID == pid) |> slice(1)
      idr_info <- idr_expanded |> filter(Protein.ID == pid) |>
        arrange(logFC) |> slice(1)  # Get the one closest to 0 but negative

      results <- bind_rows(results, tibble(
        Protein.ID = pid,
        Gene = struct_info$Gene,
        Struct_site = struct_info$site_number,
        Struct_logFC = struct_info$logFC,
        Struct_pLDDT = struct_info$pLDDT,
        IDR_site = idr_info$site_number,
        IDR_logFC = idr_info$logFC,
        IDR_pLDDT = idr_info$pLDDT
      ))
    }

    results <- results |> arrange(desc(Struct_logFC))
    print(results, n = 20)
    cat("\n")
  } else {
    cat("No matching proteins\n\n")
  }
}
