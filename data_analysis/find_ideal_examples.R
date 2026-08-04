# Find ideal examples for Figure 6E/6F
library(tidyverse)
source('data_source.R')

# Load pre-computed features
features_analysis <- read_csv(
  paste0(source_file_path, 'site_features/OGlcNAc_site_features.csv'),
  show_col_types = FALSE
)

# Filter for HEK293T only
hek_features <- features_analysis |>
  filter(cell == "HEK293T", !is.na(is_IDR), !is.na(logFC))

cat("Total HEK293T sites with IDR info:", nrow(hek_features), "\n\n")

# Distribution check
cat("=== DISTRIBUTION OF SITES ===\n")
cat("Structured sites:", sum(hek_features$is_IDR == 0), "\n")
cat("IDR sites:", sum(hek_features$is_IDR == 1), "\n\n")

cat("Structured logFC range:", range(hek_features$logFC[hek_features$is_IDR == 0]), "\n")
cat("IDR logFC range:", range(hek_features$logFC[hek_features$is_IDR == 1]), "\n\n")

# Find proteins with sites in BOTH regions
proteins_with_both <- hek_features |>
  group_by(Protein.ID, Gene) |>
  summarize(
    n_structured = sum(is_IDR == 0),
    n_IDR = sum(is_IDR == 1),
    n_total = n(),
    mean_logFC_structured = mean(logFC[is_IDR == 0], na.rm = TRUE),
    mean_logFC_IDR = mean(logFC[is_IDR == 1], na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(n_structured >= 1, n_IDR >= 1)

cat("=== ALL PROTEINS WITH SITES IN BOTH REGIONS (HEK293T) ===\n")
cat("Total proteins:", nrow(proteins_with_both), "\n\n")

print(proteins_with_both, n = Inf)

# Show detailed site data for each
cat("\n\n=== DETAILED SITE-LEVEL DATA ===\n")

for (i in 1:nrow(proteins_with_both)) {
  pid <- proteins_with_both$Protein.ID[i]
  gene <- proteins_with_both$Gene[i]

  cat("\n")
  cat(rep("=", 50), "\n", sep = "")
  cat("Protein:", pid, "(", gene, ")\n")
  cat(rep("=", 50), "\n", sep = "")

  site_details <- hek_features |>
    filter(Protein.ID == pid) |>
    dplyr::select(site_number, logFC, is_IDR, pLDDT, secondary_structure_simple) |>
    mutate(
      region = ifelse(is_IDR == 1, "IDR", "Structured"),
      site_label = paste0("S", site_number)
    ) |>
    arrange(is_IDR, desc(logFC))

  print(site_details, n = Inf)

  cat("\nSummary: Structured mean logFC =",
      round(mean(site_details$logFC[site_details$is_IDR == 0]), 3),
      ", IDR mean logFC =",
      round(mean(site_details$logFC[site_details$is_IDR == 1]), 3), "\n")
}

# ============================================
# CHECK ALL CELLS - Maybe other cells have better examples
# ============================================

cat("\n\n")
cat(rep("=", 60), "\n", sep = "")
cat("=== CHECKING ALL CELL TYPES ===\n")
cat(rep("=", 60), "\n", sep = "")

all_features <- features_analysis |>
  filter(!is.na(is_IDR), !is.na(logFC))

all_proteins_both <- all_features |>
  group_by(Protein.ID, Gene, cell) |>
  summarize(
    n_structured = sum(is_IDR == 0),
    n_IDR = sum(is_IDR == 1),
    n_total = n(),
    mean_logFC_structured = mean(logFC[is_IDR == 0], na.rm = TRUE),
    mean_logFC_IDR = mean(logFC[is_IDR == 1], na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(n_structured >= 1, n_IDR >= 1) |>
  mutate(diff = mean_logFC_structured - mean_logFC_IDR)

cat("\nProteins with sites in both regions by cell type:\n")
print(table(all_proteins_both$cell))

cat("\n=== ALL CANDIDATES ACROSS ALL CELLS ===\n")
all_proteins_both |>
  arrange(desc(diff)) |>
  print(n = 30)

# Find ideal examples across all cells
cat("\n=== IDEAL EXAMPLES (Structured ~1, IDR ~0) - ALL CELLS ===\n")
ideal_all <- all_proteins_both |>
  filter(
    mean_logFC_structured >= 0.5 & mean_logFC_structured <= 1.5,
    mean_logFC_IDR >= -0.5 & mean_logFC_IDR <= 0.3
  ) |>
  arrange(desc(diff))

if (nrow(ideal_all) > 0) {
  print(ideal_all, n = 20)
} else {
  cat("No ideal examples found.\n")
  cat("\nClosest matches (Structured > 0.3, IDR < 0.3):\n")
  close_matches <- all_proteins_both |>
    filter(mean_logFC_structured > 0.3, mean_logFC_IDR < 0.3) |>
    arrange(desc(diff))
  print(close_matches, n = 20)
}
