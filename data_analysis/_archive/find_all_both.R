# Find ALL proteins with sites in BOTH structured and IDR regions
library(tidyverse)
source('data_source.R')

features_analysis <- read_csv(
  paste0(source_file_path, 'site_features/OGlcNAc_site_features.csv'),
  show_col_types = FALSE
)

cat("=== ALL PROTEINS WITH SITES IN BOTH REGIONS ===\n\n")

for (cell_type in c("HEK293T", "HepG2", "Jurkat")) {

  cat(rep("=", 60), "\n", sep = "")
  cat("Cell Type:", cell_type, "\n")
  cat(rep("=", 60), "\n\n")

  cell_features <- features_analysis |>
    filter(cell == cell_type, !is.na(is_IDR), !is.na(logFC))

  # Find proteins with both
  proteins_both <- cell_features |>
    group_by(Protein.ID, Gene) |>
    summarize(
      n_structured = sum(is_IDR == 0),
      n_IDR = sum(is_IDR == 1),
      .groups = "drop"
    ) |>
    filter(n_structured >= 1, n_IDR >= 1)

  cat("Total proteins with both:", nrow(proteins_both), "\n\n")

  if (nrow(proteins_both) > 0) {
    for (i in 1:nrow(proteins_both)) {
      pid <- proteins_both$Protein.ID[i]
      gene <- proteins_both$Gene[i]

      cat(">>> ", pid, " (", gene, ") <<<\n", sep = "")

      sites <- cell_features |>
        filter(Protein.ID == pid) |>
        select(site_number, logFC, pLDDT, is_IDR) |>
        mutate(region = ifelse(is_IDR == 1, "IDR", "Structured")) |>
        arrange(is_IDR, desc(logFC))

      print(sites, n = Inf)
      cat("\n")
    }
  }
}

# ============================================
# Check if any protein has structured site (positive) + IDR site (any)
# ============================================

cat("\n")
cat(rep("=", 60), "\n", sep = "")
cat("=== PROTEINS: Structured (logFC > 0) + ANY IDR site ===\n")
cat(rep("=", 60), "\n\n")

all_features <- features_analysis |>
  filter(!is.na(is_IDR), !is.na(logFC))

# All structured sites with positive logFC
struct_positive <- all_features |>
  filter(is_IDR == 0, logFC > 0) |>
  distinct(Protein.ID)

# All proteins with IDR sites
has_idr <- all_features |>
  filter(is_IDR == 1) |>
  distinct(Protein.ID)

# Intersection
both <- intersect(struct_positive$Protein.ID, has_idr$Protein.ID)

cat("Proteins with structured (logFC>0) AND any IDR site:", length(both), "\n\n")

for (pid in both) {
  sites <- all_features |>
    filter(Protein.ID == pid) |>
    select(Protein.ID, Gene, cell, site_number, logFC, pLDDT, is_IDR) |>
    mutate(region = ifelse(is_IDR == 1, "IDR", "Structured")) |>
    arrange(cell, is_IDR, desc(logFC))

  gene <- sites$Gene[1]
  cat("\n>>> ", pid, " (", gene, ") <<<\n", sep = "")
  print(sites |> select(-Protein.ID, -Gene), n = Inf)
}
