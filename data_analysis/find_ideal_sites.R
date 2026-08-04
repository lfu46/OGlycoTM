# Find ideal site-level examples for Figure 6E/6F
# Criteria (individual sites, not averages):
# - At least 1 structured site with logFC ~+1 (0.7 to 1.3)
# - At least 1 IDR site with logFC slightly < 0 (-0.3 to 0)

library(tidyverse)
source('data_source.R')

# Load pre-computed features
features_analysis <- read_csv(
  paste0(source_file_path, 'site_features/OGlcNAc_site_features.csv'),
  show_col_types = FALSE
)

cat("=== SITE-LEVEL SEARCH ===\n\n")

# Check all cells
for (cell_type in c("HEK293T", "HepG2", "Jurkat")) {

  cat(rep("=", 60), "\n", sep = "")
  cat("Cell Type:", cell_type, "\n")
  cat(rep("=", 60), "\n\n")

  cell_features <- features_analysis |>
    filter(cell == cell_type, !is.na(is_IDR), !is.na(logFC))

  # Find structured sites with logFC ~1 (0.7 to 1.3)
  structured_ideal <- cell_features |>
    filter(is_IDR == 0, logFC >= 0.7, logFC <= 1.3) |>
    select(Protein.ID, Gene, site_number, logFC, pLDDT, is_IDR) |>
    mutate(region = "Structured")

  cat("Structured sites with logFC 0.7-1.3:\n")
  if (nrow(structured_ideal) > 0) {
    print(structured_ideal, n = Inf)
  } else {
    cat("None found\n")
  }

  # Find IDR sites with logFC slightly below 0 (-0.3 to 0)
  idr_ideal <- cell_features |>
    filter(is_IDR == 1, logFC >= -0.3, logFC <= 0) |>
    select(Protein.ID, Gene, site_number, logFC, pLDDT, is_IDR) |>
    mutate(region = "IDR")

  cat("\nIDR sites with logFC -0.3 to 0:\n")
  if (nrow(idr_ideal) > 0) {
    print(idr_ideal, n = Inf)
  } else {
    cat("None found\n")
  }

  # Find proteins that have BOTH types of sites
  if (nrow(structured_ideal) > 0 && nrow(idr_ideal) > 0) {
    common_proteins <- intersect(structured_ideal$Protein.ID, idr_ideal$Protein.ID)

    cat("\n*** PROTEINS WITH BOTH IDEAL SITE TYPES ***\n")
    if (length(common_proteins) > 0) {
      for (pid in common_proteins) {
        gene <- cell_features |> filter(Protein.ID == pid) |> pull(Gene) |> unique()
        cat("\n>>> ", pid, " (", gene, ") <<<\n", sep = "")

        # Show all sites for this protein
        all_sites <- cell_features |>
          filter(Protein.ID == pid) |>
          select(site_number, logFC, pLDDT, is_IDR, secondary_structure_simple) |>
          mutate(region = ifelse(is_IDR == 1, "IDR", "Structured")) |>
          arrange(is_IDR, site_number)
        print(all_sites, n = Inf)

        # Mark which sites meet criteria
        cat("\nIdeal structured site(s): ")
        cat(structured_ideal |> filter(Protein.ID == pid) |> pull(site_number), "\n")
        cat("Ideal IDR site(s): ")
        cat(idr_ideal |> filter(Protein.ID == pid) |> pull(site_number), "\n")
      }
    } else {
      cat("No proteins found with both types.\n")
    }
  }

  cat("\n\n")
}

# ============================================
# RELAXED SEARCH: Expand criteria
# ============================================

cat(rep("=", 60), "\n", sep = "")
cat("=== RELAXED CRITERIA SEARCH ===\n")
cat("Structured: logFC 0.5 to 1.5\n")
cat("IDR: logFC -0.5 to 0.1\n")
cat(rep("=", 60), "\n\n")

for (cell_type in c("HEK293T", "HepG2", "Jurkat")) {

  cat("--- ", cell_type, " ---\n", sep = "")

  cell_features <- features_analysis |>
    filter(cell == cell_type, !is.na(is_IDR), !is.na(logFC))

  # Relaxed structured sites
  structured_relaxed <- cell_features |>
    filter(is_IDR == 0, logFC >= 0.5, logFC <= 1.5) |>
    select(Protein.ID, Gene, site_number, logFC, pLDDT)

  # Relaxed IDR sites
  idr_relaxed <- cell_features |>
    filter(is_IDR == 1, logFC >= -0.5, logFC <= 0.1) |>
    select(Protein.ID, Gene, site_number, logFC, pLDDT)

  # Find common proteins (excluding TAB2)
  common_proteins <- setdiff(
    intersect(structured_relaxed$Protein.ID, idr_relaxed$Protein.ID),
    "Q9NYJ8"  # Exclude TAB2
  )

  if (length(common_proteins) > 0) {
    cat("Found proteins (excluding TAB2):\n")
    for (pid in common_proteins) {
      gene <- cell_features |> filter(Protein.ID == pid) |> pull(Gene) |> unique()
      cat("\n>>> ", pid, " (", gene, ") <<<\n", sep = "")

      all_sites <- cell_features |>
        filter(Protein.ID == pid) |>
        select(site_number, logFC, pLDDT, is_IDR) |>
        mutate(region = ifelse(is_IDR == 1, "IDR", "Structured")) |>
        arrange(is_IDR, site_number)
      print(all_sites, n = Inf)

      cat("Matching structured site(s): S")
      cat(paste(structured_relaxed |> filter(Protein.ID == pid) |> pull(site_number), collapse = ", S"), "\n")
      cat("Matching IDR site(s): S")
      cat(paste(idr_relaxed |> filter(Protein.ID == pid) |> pull(site_number), collapse = ", S"), "\n")
    }
  } else {
    cat("No matching proteins (excluding TAB2)\n")
  }
  cat("\n")
}

# ============================================
# SHOW ALL AVAILABLE SITES FOR REFERENCE
# ============================================

cat("\n")
cat(rep("=", 60), "\n", sep = "")
cat("=== ALL STRUCTURED SITES (for reference) ===\n")
cat(rep("=", 60), "\n\n")

all_structured <- features_analysis |>
  filter(!is.na(is_IDR), !is.na(logFC), is_IDR == 0) |>
  select(Protein.ID, Gene, cell, site_number, logFC, pLDDT) |>
  arrange(desc(logFC))

print(all_structured, n = 50)

cat("\n\n=== IDR SITES WITH logFC NEAR 0 (-0.5 to 0.5) ===\n\n")

idr_near_zero <- features_analysis |>
  filter(!is.na(is_IDR), !is.na(logFC), is_IDR == 1, logFC >= -0.5, logFC <= 0.5) |>
  select(Protein.ID, Gene, cell, site_number, logFC, pLDDT) |>
  arrange(logFC)

print(idr_near_zero, n = 50)
