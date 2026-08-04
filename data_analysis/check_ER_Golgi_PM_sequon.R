# Check N-glycosylation sequon in O-GlcNAc PSMs from ER/Golgi/PM-only proteins
# All three cell lines: HEK293T, HepG2, Jurkat

library(tidyverse)

source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'

# Load O-GlcNAc filtered data (Level1/1b, no Cys artifacts)
OGlcNAc_HEK293T <- read_csv(paste0(source_file_path, 'filtered/OGlcNAc_HEK293T.csv'), show_col_types = FALSE)
OGlcNAc_HepG2 <- read_csv(paste0(source_file_path, 'filtered/OGlcNAc_HepG2.csv'), show_col_types = FALSE)
OGlcNAc_Jurkat <- read_csv(paste0(source_file_path, 'filtered/OGlcNAc_Jurkat.csv'), show_col_types = FALSE)

# Load HPA subcellular location
subcellular_location <- read_tsv(
  paste0(source_file_path, 'reference/subcellular_location.tsv'),
  show_col_types = FALSE
)

# =============================================================================
# Filter for proteins with ER/Golgi/PM as ONLY location
# =============================================================================
ER_Golgi_PM_targets <- c("Endoplasmic reticulum", "Golgi apparatus", "Plasma membrane")

hpa_erpm_only <- subcellular_location |>
  filter(!str_detect(`Main location`, ";")) |>
  filter(Reliability != "Uncertain") |>
  filter(`Main location` %in% ER_Golgi_PM_targets) |>
  filter(is.na(`Additional location`) | `Additional location` == "") |>
  dplyr::select(Gene_name = `Gene name`, Location = `Main location`, Reliability) |>
  distinct(Gene_name, Location, .keep_all = TRUE)

# Remove genes with multiple entries mapping to different locations
genes_multiple <- hpa_erpm_only |>
  group_by(Gene_name) |> filter(n() > 1) |> pull(Gene_name) |> unique()
hpa_erpm_only <- hpa_erpm_only |> filter(!(Gene_name %in% genes_multiple))

cat("=== HPA proteins with ER/Golgi/PM as ONLY location ===\n")
print(table(hpa_erpm_only$Location))
cat("Total:", nrow(hpa_erpm_only), "\n\n")

# =============================================================================
# Function to analyze one cell line
# =============================================================================
analyze_cell <- function(df, cell_name, hpa) {
  cat("===============================================================\n")
  cat("  ", cell_name, "\n")
  cat("===============================================================\n")

  merged <- df |> inner_join(hpa, by = c("Gene" = "Gene_name"))
  cat("Total O-GlcNAc PSMs in ER/Golgi/PM-only proteins:", nrow(merged), "\n")

  localized <- merged |> filter(Confidence.Level %in% c("Level1", "Level1b"))
  cat("Localized PSMs (Level1/1b):", nrow(localized), "\n")
  cat("Unique proteins:", n_distinct(localized$Gene), "\n")

  if (nrow(localized) == 0) {
    cat("No localized PSMs found.\n\n")
    return(NULL)
  }

  cat("\nBy location:\n")
  loc_summary <- localized |> group_by(Location) |> summarise(
    n_PSMs = n(), n_proteins = n_distinct(Gene),
    n_sequon = sum(Has.N.Glyc.Sequon == TRUE, na.rm = TRUE),
    pct_sequon = round(100 * n_sequon / n(), 1), .groups = "drop"
  )
  print(loc_summary)

  cat("\nN-glyc sequon breakdown:\n")
  print(table(localized$Has.N.Glyc.Sequon, useNA = "always"))

  # Unique sites
  sites <- localized |>
    mutate(
      site_residue = str_extract(Site.Probabilities, "\\d+"),
      site_index = paste0(Protein.ID, "_", site_residue)
    ) |>
    group_by(Gene, Protein.ID, Location, site_index) |>
    summarise(
      n_PSMs = n(),
      has_sequon = any(Has.N.Glyc.Sequon == TRUE),
      best_confidence = ifelse(any(Confidence.Level == "Level1"), "Level1", "Level1b"),
      example_peptide = first(Peptide),
      example_mods = first(Assigned.Modifications),
      example_site_prob = first(Site.Probabilities),
      .groups = "drop"
    )

  cat("\nUnique sites:", nrow(sites), "\n")
  cat("Sites with N-glyc sequon:", sum(sites$has_sequon), "\n")
  cat("Sites without N-glyc sequon:", sum(!sites$has_sequon), "\n\n")
  print(sites |> arrange(Location, Gene), n = Inf, width = 200)

  # Comparison
  overall_sequon <- sum(df$Has.N.Glyc.Sequon == TRUE, na.rm = TRUE)
  overall_total <- nrow(df)
  erpm_sequon <- sum(localized$Has.N.Glyc.Sequon == TRUE, na.rm = TRUE)
  erpm_total <- nrow(localized)
  cat("\nSequon rate - All O-GlcNAc:", overall_sequon, "/", overall_total,
      "(", round(100 * overall_sequon / overall_total, 1), "%)\n")
  cat("Sequon rate - ER/Golgi/PM-only:", erpm_sequon, "/", erpm_total,
      "(", round(100 * erpm_sequon / erpm_total, 1), "%)\n\n")

  return(sites)
}

# =============================================================================
# Run for all three cell lines
# =============================================================================
sites_HEK293T <- analyze_cell(OGlcNAc_HEK293T, "HEK293T", hpa_erpm_only)
sites_HepG2   <- analyze_cell(OGlcNAc_HepG2, "HepG2", hpa_erpm_only)
sites_Jurkat  <- analyze_cell(OGlcNAc_Jurkat, "Jurkat", hpa_erpm_only)

# =============================================================================
# Combined summary
# =============================================================================
cat("\n===============================================================\n")
cat("  COMBINED SUMMARY\n")
cat("===============================================================\n")

all_sites <- bind_rows(
  if (!is.null(sites_HEK293T)) sites_HEK293T |> mutate(cell = "HEK293T"),
  if (!is.null(sites_HepG2))   sites_HepG2   |> mutate(cell = "HepG2"),
  if (!is.null(sites_Jurkat))  sites_Jurkat  |> mutate(cell = "Jurkat")
)

if (nrow(all_sites) > 0) {
  cat("\nAll unique sites across cell lines:\n")
  print(all_sites |> arrange(Location, Gene, cell), n = Inf, width = 200)

  cat("\nTotal unique sites:", nrow(all_sites), "\n")
  cat("With N-glyc sequon:", sum(all_sites$has_sequon), "\n")
  cat("Without N-glyc sequon:", sum(!all_sites$has_sequon), "\n")

  # Unique proteins across all
  cat("\nUnique proteins across all cell lines:", n_distinct(all_sites$Gene), "\n")

  # Check overlap
  cat("\nSites by cell line:\n")
  print(all_sites |> group_by(cell) |> summarise(
    n_sites = n(), n_proteins = n_distinct(Gene),
    n_sequon = sum(has_sequon), .groups = "drop"
  ))
}
