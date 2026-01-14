# Figure 5: Subcellular Localization Analysis of O-GlcNAc Proteins

library(tidyverse)

# Source data paths, colors, and differential analysis results
source('data_source_DE.R')

# =============================================================================
# Data Import - Subcellular Location Annotation
# =============================================================================
# Data source: Human Protein Atlas subcellular location

# Load subcellular location data
subcellular_location <- read_tsv(
  paste0(source_file_path, 'reference/subcellular_location.tsv')
)

# Filter for:
# 1. Single location in Main location (no ";" separator)
# 2. Reliability not "Uncertain"
subcellular_location_filtered <- subcellular_location |>
  filter(!str_detect(`Main location`, ";")) |>
  filter(Reliability != "Uncertain") |>
  dplyr::select(Gene_name = `Gene name`, Location = `Main location`, Reliability)

cat("Subcellular location data after filtering:\n")
cat("Total proteins with single location and reliable annotation:", nrow(subcellular_location_filtered), "\n")

# Check if each gene has only one distinct location annotation
duplicate_genes <- subcellular_location_filtered |>
  group_by(Gene_name) |>
  summarise(
    n_entries = n(),
    n_locations = n_distinct(Location),
    locations = paste(unique(Location), collapse = "; "),
    .groups = "drop"
  ) |>
  filter(n_entries > 1 | n_locations > 1)

if (nrow(duplicate_genes) > 0) {
  cat("\nWARNING: Found genes with multiple entries or locations:\n")
  print(duplicate_genes, n = Inf)

  # Keep only distinct gene-location pairs (remove exact duplicates)
  subcellular_location_filtered <- subcellular_location_filtered |>
    distinct(Gene_name, Location, .keep_all = TRUE)

  # Check again for genes with multiple different locations
  genes_multiple_locations <- subcellular_location_filtered |>
    group_by(Gene_name) |>
    filter(n() > 1) |>
    ungroup()

  if (nrow(genes_multiple_locations) > 0) {
    cat("\nGenes with multiple different locations (will be removed):\n")
    print(genes_multiple_locations, n = Inf)

    # Remove genes with multiple different locations
    genes_to_remove <- genes_multiple_locations |> pull(Gene_name) |> unique()
    subcellular_location_filtered <- subcellular_location_filtered |>
      filter(!(Gene_name %in% genes_to_remove))

    cat("\nRemoved", length(genes_to_remove), "genes with ambiguous locations\n")
  }
} else {
  cat("\nAll genes have unique location annotations.\n")
}

cat("Final filtered proteins:", nrow(subcellular_location_filtered), "\n")
cat("Unique genes:", n_distinct(subcellular_location_filtered$Gene_name), "\n")

# Create output directory for intermediate files
dir.create(paste0(source_file_path, 'subcellular_location'), showWarnings = FALSE)

# Save filtered subcellular location data
write_csv(
  subcellular_location_filtered,
  paste0(source_file_path, 'subcellular_location/subcellular_location_filtered.csv')
)
cat("\nSaved: subcellular_location_filtered.csv\n")

# Check location distribution
cat("\nLocation distribution:\n")
subcellular_location_filtered |>
  group_by(Location) |>
  summarise(n = n(), .groups = "drop") |>
  arrange(desc(n)) |>
  print(n = 30)

# =============================================================================
# Figure 5A - Proportion Barplot of Subcellular Locations
# =============================================================================
# Stacked barplot showing proportion of O-GlcNAc proteins in each location
# Locations with <10 proteins or no annotation go to "Other"

# Load required libraries
library(tidyverse)

# Source data paths, colors, and differential analysis results
source('data_source_DE.R')

# Create output directory
dir.create(paste0(source_file_path, 'subcellular_location'), showWarnings = FALSE)

# Load subcellular location data
subcellular_location_filtered <- read_csv(
  paste0(source_file_path, 'subcellular_location/subcellular_location_filtered.csv')
)
cat("Loaded: subcellular_location_filtered.csv\n")

# Function to prepare data for proportion barplot
prepare_proportion_data <- function(de_data, location_data, cell_name, min_n = 10) {
  # Left join to keep all proteins (including those without annotation)
  data_with_loc <- de_data |>
    left_join(location_data, by = c("Gene" = "Gene_name")) |>
    dplyr::select(Protein.ID, Gene, Location)

  # Count proteins per location

  location_counts <- data_with_loc |>
    mutate(Location = ifelse(is.na(Location), "No annotation", Location)) |>
    group_by(Location) |>
    summarise(n = n(), .groups = "drop")

  # Identify locations with <10 proteins
  small_locations <- location_counts |>
    filter(n < min_n) |>
    pull(Location)

  # Reclassify small groups and no annotation to "Other"
  data_final <- data_with_loc |>
    mutate(
      Location = case_when(
        is.na(Location) ~ "Other",
        Location %in% small_locations ~ "Other",
        TRUE ~ Location
      )
    ) |>
    group_by(Location) |>
    summarise(n = n(), .groups = "drop") |>
    mutate(
      cell = cell_name,
      proportion = n / sum(n)
    )

  return(data_final)
}

# Prepare data for each cell type
prop_HEK293T <- prepare_proportion_data(OGlcNAc_protein_DE_HEK293T, subcellular_location_filtered, "HEK293T")
prop_HepG2 <- prepare_proportion_data(OGlcNAc_protein_DE_HepG2, subcellular_location_filtered, "HepG2")
prop_Jurkat <- prepare_proportion_data(OGlcNAc_protein_DE_Jurkat, subcellular_location_filtered, "Jurkat")

# Combine data
prop_combined <- bind_rows(prop_HEK293T, prop_HepG2, prop_Jurkat) |>
  mutate(cell = factor(cell, levels = c("HEK293T", "HepG2", "Jurkat")))

# Check the data
cat("\nProportion data for barplot:\n")
prop_combined |>
  pivot_wider(names_from = cell, values_from = c(n, proportion), values_fill = 0) |>
  print(n = Inf)

# Save proportion data
write_csv(
  prop_combined,
  paste0(source_file_path, 'subcellular_location/OGlcNAc_protein_location_proportion.csv')
)
cat("Saved: OGlcNAc_protein_location_proportion.csv\n")

# Order locations by total count (descending), with "Other" at the end
location_order <- prop_combined |>
  group_by(Location) |>
  summarise(total = sum(n), .groups = "drop") |>
  filter(Location != "Other") |>
  arrange(desc(total)) |>
  pull(Location)
location_order <- c(location_order, "Other")

prop_combined <- prop_combined |>
  mutate(Location = factor(Location, levels = rev(location_order)))

# Define publication-quality color palette for subcellular locations
# Compatible with colors from data_source.R, with grey for "Other"
colors_location <- c(
  "Nucleoplasm" = "#3C5488",
  "Cytosol" = "#00A087",
  "Nucleoli" = "#4DBBD5",
  "Mitochondria" = "#E64B35",
  "Vesicles" = "#F39B7F",
  "Golgi apparatus" = "#8491B4",
  "Endoplasmic reticulum" = "#7E6148",
  "Nuclear membrane" = "#B09C85",
  "Plasma membrane" = "#91D1C2",
  "Nuclear bodies" = "#DC0000",
  "Nuclear speckles" = "#7AA6DC",
  "Centrosome" = "#868686",
  "Other" = "#BEBEBE"
)

# Create proportion barplot
figure5A <- prop_combined |>
  ggplot(aes(x = cell, y = proportion, fill = Location)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = colors_location, na.value = "#BEBEBE") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  labs(x = "", y = "Proportion of O-GlcNAc proteins", fill = "Location") +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(color = "black", size = 8, angle = 30, hjust = 1),
    axis.text.y = element_text(color = "black", size = 8),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.2, "cm"),
    legend.spacing.y = unit(0.05, "cm")
  )

# Create output directory
dir.create(paste0(figure_file_path, "Figure5"), showWarnings = FALSE)

ggsave(
  filename = paste0(figure_file_path, 'Figure5/Figure5A.pdf'),
  plot = figure5A,
  height = 1.8, width = 2.5, units = 'in'
)

cat("\nFigure 5A proportion barplot saved to:", figure_file_path, "Figure5/Figure5A.pdf\n")

# =============================================================================
# Figure 5B - Wilcoxon Tests and Violin Boxplot
# =============================================================================
# Compare logFC distributions across subcellular locations and cell types

# Load required libraries
library(tidyverse)
library(ggpubr)

# Source data paths, colors, and differential analysis results
source('data_source_DE.R')

# Create output directories
dir.create(paste0(source_file_path, 'subcellular_location'), showWarnings = FALSE)
dir.create(paste0(figure_file_path, "Figure5"), showWarnings = FALSE)

# Load subcellular location data
subcellular_location_filtered <- read_csv(
  paste0(source_file_path, 'subcellular_location/subcellular_location_filtered.csv')
)
cat("Loaded: subcellular_location_filtered.csv\n")

# Merge with DE data for each cell type
merge_with_location <- function(de_data, cell_name) {
  de_data |>
    inner_join(subcellular_location_filtered, by = c("Gene" = "Gene_name")) |>
    dplyr::select(Protein.ID, Gene, Location, logFC, adj.P.Val) |>
    mutate(cell = cell_name)
}

OGlcNAc_loc_HEK293T <- merge_with_location(OGlcNAc_protein_DE_HEK293T, "HEK293T")
OGlcNAc_loc_HepG2 <- merge_with_location(OGlcNAc_protein_DE_HepG2, "HepG2")
OGlcNAc_loc_Jurkat <- merge_with_location(OGlcNAc_protein_DE_Jurkat, "Jurkat")

cat("\nO-GlcNAc proteins with location annotation:\n")
cat("HEK293T:", nrow(OGlcNAc_loc_HEK293T), "\n")
cat("HepG2:", nrow(OGlcNAc_loc_HepG2), "\n")
cat("Jurkat:", nrow(OGlcNAc_loc_Jurkat), "\n")

# Combine all data
OGlcNAc_location_combined <- bind_rows(
  OGlcNAc_loc_HEK293T,
  OGlcNAc_loc_HepG2,
  OGlcNAc_loc_Jurkat
)

# Save combined location data
write_csv(
  OGlcNAc_location_combined,
  paste0(source_file_path, 'subcellular_location/OGlcNAc_protein_location_combined.csv')
)
cat("\nSaved: OGlcNAc_protein_location_combined.csv\n")

# Check protein counts per location per cell type
cat("\nProtein counts per location per cell type:\n")
location_counts <- OGlcNAc_location_combined |>
  group_by(cell, Location) |>
  summarise(n = n(), .groups = "drop") |>
  pivot_wider(names_from = cell, values_from = n, values_fill = 0)
print(location_counts, n = 30)

# Save location counts
write_csv(
  location_counts,
  paste0(source_file_path, 'subcellular_location/OGlcNAc_protein_location_counts.csv')
)
cat("Saved: OGlcNAc_protein_location_counts.csv\n")

# Filter locations with at least 10 proteins in at least one cell type
locations_with_min_10 <- OGlcNAc_location_combined |>
  group_by(cell, Location) |>
  summarise(n = n(), .groups = "drop") |>
  filter(n >= 10) |>
  pull(Location) |>
  unique()

cat("\nLocations with >= 10 proteins in at least one cell type:\n")
print(locations_with_min_10)

# -----------------------------------------------------------------------------
# Wilcoxon Test 1: Across Locations Within Each Cell Type
# -----------------------------------------------------------------------------
cat("\n\n========== Wilcoxon Test: Across Locations Within Each Cell Type ==========\n")

wilcox_across_locations <- function(data, cell_name, min_n = 10) {
  # Filter data for this cell type
  cell_data <- data |> filter(cell == cell_name)

  # Get locations with at least min_n proteins
  valid_locations <- cell_data |>
    group_by(Location) |>
    summarise(n = n(), .groups = "drop") |>
    filter(n >= min_n) |>
    pull(Location)

  if (length(valid_locations) < 2) {
    return(tibble())
  }

  # Generate all pairs of locations
  location_pairs <- combn(valid_locations, 2, simplify = FALSE)

  # Perform Wilcoxon test for each pair
  results <- map_dfr(location_pairs, function(pair) {
    loc1_data <- cell_data |> filter(Location == pair[1]) |> pull(logFC)
    loc2_data <- cell_data |> filter(Location == pair[2]) |> pull(logFC)

    wilcox_result <- wilcox.test(loc1_data, loc2_data)

    tibble(
      cell = cell_name,
      Location_1 = pair[1],
      Location_2 = pair[2],
      n_Location_1 = length(loc1_data),
      n_Location_2 = length(loc2_data),
      mean_logFC_1 = mean(loc1_data, na.rm = TRUE),
      mean_logFC_2 = mean(loc2_data, na.rm = TRUE),
      W_statistic = wilcox_result$statistic,
      p_value = wilcox_result$p.value
    )
  })

  return(results)
}

# Run Wilcoxon tests for each cell type
wilcox_locations_HEK293T <- wilcox_across_locations(OGlcNAc_location_combined, "HEK293T")
wilcox_locations_HepG2 <- wilcox_across_locations(OGlcNAc_location_combined, "HepG2")
wilcox_locations_Jurkat <- wilcox_across_locations(OGlcNAc_location_combined, "Jurkat")

# Combine results
wilcox_across_locations_results <- bind_rows(
  wilcox_locations_HEK293T,
  wilcox_locations_HepG2,
  wilcox_locations_Jurkat
) |>
  mutate(p_adj = p.adjust(p_value, method = "BH")) |>
  arrange(p_value)

cat("\nWilcoxon Test Results - Across Locations Within Each Cell Type:\n")
print(wilcox_across_locations_results, n = Inf)

# Save results
write_csv(
  wilcox_across_locations_results,
  paste0(source_file_path, 'subcellular_location/Wilcoxon_across_locations_within_cell.csv')
)

# -----------------------------------------------------------------------------
# Wilcoxon Test 2: Across Cell Types Within Each Location
# -----------------------------------------------------------------------------
cat("\n\n========== Wilcoxon Test: Across Cell Types Within Each Location ==========\n")

wilcox_across_cells <- function(data, location_name, min_n = 10) {
  # Filter data for this location
  loc_data <- data |> filter(Location == location_name)

  # Get cell types with at least min_n proteins
  valid_cells <- loc_data |>
    group_by(cell) |>
    summarise(n = n(), .groups = "drop") |>
    filter(n >= min_n) |>
    pull(cell)

  if (length(valid_cells) < 2) {
    return(tibble())
  }

  # Generate all pairs of cell types
  cell_pairs <- combn(valid_cells, 2, simplify = FALSE)

  # Perform Wilcoxon test for each pair
  results <- map_dfr(cell_pairs, function(pair) {
    cell1_data <- loc_data |> filter(cell == pair[1]) |> pull(logFC)
    cell2_data <- loc_data |> filter(cell == pair[2]) |> pull(logFC)

    wilcox_result <- wilcox.test(cell1_data, cell2_data)

    tibble(
      Location = location_name,
      Cell_1 = pair[1],
      Cell_2 = pair[2],
      n_Cell_1 = length(cell1_data),
      n_Cell_2 = length(cell2_data),
      mean_logFC_1 = mean(cell1_data, na.rm = TRUE),
      mean_logFC_2 = mean(cell2_data, na.rm = TRUE),
      W_statistic = wilcox_result$statistic,
      p_value = wilcox_result$p.value
    )
  })

  return(results)
}

# Get all unique locations
all_locations <- unique(OGlcNAc_location_combined$Location)

# Run Wilcoxon tests for each location
wilcox_across_cells_results <- map_dfr(all_locations, function(loc) {
  wilcox_across_cells(OGlcNAc_location_combined, loc)
}) |>
  mutate(p_adj = p.adjust(p_value, method = "BH")) |>
  arrange(p_value)

cat("\nWilcoxon Test Results - Across Cell Types Within Each Location:\n")
print(wilcox_across_cells_results, n = Inf)

# Save results
write_csv(
  wilcox_across_cells_results,
  paste0(source_file_path, 'subcellular_location/Wilcoxon_across_cells_within_location.csv')
)

# -----------------------------------------------------------------------------
# Summary Statistics
# -----------------------------------------------------------------------------
cat("\n\n========== Summary Statistics ==========\n")

# Summary: Mean logFC per location per cell type
cat("\nMean logFC per Location per Cell Type:\n")
summary_logFC <- OGlcNAc_location_combined |>
  group_by(Location, cell) |>
  summarise(
    n = n(),
    mean_logFC = mean(logFC, na.rm = TRUE),
    sd_logFC = sd(logFC, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(n >= 10) |>
  arrange(Location, cell)
print(summary_logFC, n = Inf)

# Save summary
write_csv(
  summary_logFC,
  paste0(source_file_path, 'subcellular_location/logFC_summary_by_location_cell.csv')
)

# Significant comparisons summary
cat("\nSignificant comparisons (p_adj < 0.05) - Across Locations:\n")
wilcox_across_locations_results |>
  filter(p_adj < 0.05) |>
  print(n = Inf)

cat("\nSignificant comparisons (p_adj < 0.05) - Across Cell Types:\n")
wilcox_across_cells_results |>
  filter(p_adj < 0.05) |>
  print(n = Inf)

# -----------------------------------------------------------------------------
# Figure 5B Dotplot - logFC by Subcellular Location
# -----------------------------------------------------------------------------

# Load required libraries
library(tidyverse)
library(ggpubr)

# Source data paths, colors, and differential analysis results
source('data_source_DE.R')

# Create output directory
dir.create(paste0(figure_file_path, "Figure5"), showWarnings = FALSE)

# Load pre-computed data
OGlcNAc_location_combined <- read_csv(
  paste0(source_file_path, 'subcellular_location/OGlcNAc_protein_location_combined.csv')
)
wilcox_across_locations_results <- read_csv(
  paste0(source_file_path, 'subcellular_location/Wilcoxon_across_locations_within_cell.csv')
)
cat("Loaded: OGlcNAc_protein_location_combined.csv and Wilcoxon results\n")

# Define color palette for subcellular locations
colors_location <- c(
  "Nucleoplasm" = "#3C5488",
  "Cytosol" = "#00A087",
  "Nucleoli" = "#4DBBD5",
  "Mitochondria" = "#E64B35",
  "Vesicles" = "#F39B7F",
  "Golgi apparatus" = "#8491B4",
  "Endoplasmic reticulum" = "#7E6148",
  "Nuclear membrane" = "#B09C85",
  "Plasma membrane" = "#91D1C2",
  "Nuclear bodies" = "#DC0000",
  "Nuclear speckles" = "#7AA6DC",
  "Centrosome" = "#868686",
  "Other" = "#BEBEBE"
)

# Filter for HEK293T and HepG2 only (Jurkat has no significant differences)
# Use locations with >= 10 proteins (already filtered in Wilcoxon test)
valid_locations_for_plot <- wilcox_across_locations_results |>
  filter(cell %in% c("HEK293T", "HepG2")) |>
  dplyr::select(Location_1, Location_2) |>
  unlist() |>
  unique()

cat("\nLocations for dotplot:", paste(valid_locations_for_plot, collapse = ", "), "\n")

# Abbreviation mapping for long location names
location_abbrev <- c(
  "Endoplasmic reticulum" = "ER",
  "Golgi apparatus" = "Golgi",
  "Nuclear speckles" = "Nuc. speckles",
  "Plasma membrane" = "PM"
)

# Prepare data for dotplot
dotplot_data <- OGlcNAc_location_combined |>
  filter(cell %in% c("HEK293T", "HepG2")) |>
  filter(Location %in% valid_locations_for_plot) |>
  mutate(
    Location_abbrev = ifelse(Location %in% names(location_abbrev),
                             location_abbrev[Location], Location),
    cell = factor(cell, levels = c("HEK293T", "HepG2")),
    Location_abbrev = factor(Location_abbrev)
  )

cat("Proteins for dotplot:", nrow(dotplot_data), "\n")

# Prepare significant comparisons for statistical annotations
sig_comparisons <- wilcox_across_locations_results |>
  filter(p_adj < 0.05) |>
  filter(cell %in% c("HEK293T", "HepG2")) |>
  dplyr::select(cell, Location_1, Location_2, p_adj)

cat("\nSignificant comparisons for annotations:\n")
print(sig_comparisons)

# Create comparison list for each cell type
comparisons_HEK293T <- sig_comparisons |>
  filter(cell == "HEK293T") |>
  rowwise() |>
  mutate(comparison = list(c(Location_1, Location_2))) |>
  pull(comparison)

comparisons_HepG2 <- sig_comparisons |>
  filter(cell == "HepG2") |>
  rowwise() |>
  mutate(comparison = list(c(Location_1, Location_2))) |>
  pull(comparison)

# Get unique locations involved in significant comparisons for each cell type
sig_locations_HEK293T <- sig_comparisons |>
  filter(cell == "HEK293T") |>
  dplyr::select(Location_1, Location_2) |>
  unlist() |>
  unique()

sig_locations_HepG2 <- sig_comparisons |>
  filter(cell == "HepG2") |>
  dplyr::select(Location_1, Location_2) |>
  unlist() |>
  unique()

cat("\nHEK293T significant locations (", length(sig_locations_HEK293T), "):", paste(sig_locations_HEK293T, collapse = ", "), "\n")
cat("HepG2 significant locations (", length(sig_locations_HepG2), "):", paste(sig_locations_HepG2, collapse = ", "), "\n")

# Prepare stat annotations for each cell type
prepare_stat_annotations <- function(sig_data, cell_name) {
  sig_data |>
    filter(cell == cell_name) |>
    mutate(
      group1 = ifelse(Location_1 %in% names(location_abbrev),
                      location_abbrev[Location_1], Location_1),
      group2 = ifelse(Location_2 %in% names(location_abbrev),
                      location_abbrev[Location_2], Location_2),
      p.signif = case_when(
        p_adj < 0.0001 ~ "****",
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) |>
    arrange(p_adj) |>
    mutate(y.position = 2.5 + (row_number() - 1) * 0.35)
}

stat_HEK293T <- prepare_stat_annotations(sig_comparisons, "HEK293T")
stat_HepG2 <- prepare_stat_annotations(sig_comparisons, "HepG2")

cat("\nHEK293T: ", nrow(stat_HEK293T), " significant pairs\n")
cat("HepG2: ", nrow(stat_HepG2), " significant pairs\n")

# Filter dotplot data to only include locations in significant comparisons
dotplot_HEK293T <- OGlcNAc_location_combined |>
  filter(cell == "HEK293T") |>
  filter(Location %in% sig_locations_HEK293T) |>
  mutate(
    Location_abbrev = ifelse(Location %in% names(location_abbrev),
                             location_abbrev[Location], Location),
    Location_abbrev = factor(Location_abbrev)
  )

dotplot_HepG2 <- OGlcNAc_location_combined |>
  filter(cell == "HepG2") |>
  filter(Location %in% sig_locations_HepG2) |>
  mutate(
    Location_abbrev = ifelse(Location %in% names(location_abbrev),
                             location_abbrev[Location], Location),
    Location_abbrev = factor(Location_abbrev)
  )

cat("\nHEK293T x-axis locations:", length(unique(dotplot_HEK293T$Location_abbrev)), "\n")
cat("HepG2 x-axis locations:", length(unique(dotplot_HepG2$Location_abbrev)), "\n")

# Create plot for HEK293T
plot_HEK293T <- dotplot_HEK293T |>
  ggplot(aes(x = Location_abbrev, y = logFC, color = Location)) +
  geom_hline(yintercept = 0, color = "grey50", linetype = "solid", linewidth = 0.5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  stat_summary(fun = median, geom = "crossbar", width = 0.7, linewidth = 0.3, color = "black") +
  scale_color_manual(values = colors_location) +
  stat_pvalue_manual(
    stat_HEK293T,
    label = "p.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    size = 3,
    bracket.size = 0.2,
    tip.length = 0.005
  ) +
  labs(x = "", y = expression(log[2](Tuni/Ctrl)), title = "HEK293T") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 9, hjust = 0.5, margin = margin(2, 0, 2, 0)),
    plot.title.position = "panel",
    axis.title.y = element_text(size = 9),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 9),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# Create plot for HepG2
plot_HepG2 <- dotplot_HepG2 |>
  ggplot(aes(x = Location_abbrev, y = logFC, color = Location)) +
  geom_hline(yintercept = 0, color = "grey50", linetype = "solid", linewidth = 0.5) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  stat_summary(fun = median, geom = "crossbar", width = 0.7, linewidth = 0.3, color = "black") +
  scale_color_manual(values = colors_location) +
  stat_pvalue_manual(
    stat_HepG2,
    label = "p.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    size = 3,
    bracket.size = 0.2,
    tip.length = 0.005
  ) +
  labs(x = "", y = "", title = "HepG2") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 9, hjust = 0.5, margin = margin(2, 0, 2, 0)),
    plot.title.position = "panel",
    axis.title.y = element_text(size = 9),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 9),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# Combine plots
figure5B <- ggarrange(plot_HEK293T, plot_HepG2, ncol = 2, nrow = 1)

ggsave(
  filename = paste0(figure_file_path, 'Figure5/Figure5B.pdf'),
  plot = figure5B,
  height = 1.8, width = 3.5, units = 'in'
)

cat("\nFigure 5B dotplot saved to:", figure_file_path, "Figure5/Figure5B.pdf\n")

cat("\nFigure 5B complete.\n")
cat("Results saved to:", source_file_path, "subcellular_location/\n")

# =============================================================================
# Figure 5C-5E - Location-Specific O-GlcNAc Protein Dot Plots
# =============================================================================
# Dot plots for selected proteins from specific subcellular locations:
# - Figure 5C: Plasma Membrane proteins (HEK293T)
# - Figure 5D: ER proteins (HepG2)
# - Figure 5E: Mitochondria proteins (HepG2)

# Load required libraries
library(tidyverse)

# Source data paths, colors, and differential analysis results
source('data_source_DE.R')

# Create output directories
dir.create(paste0(source_file_path, 'subcellular_location'), showWarnings = FALSE)
dir.create(paste0(figure_file_path, "Figure5"), showWarnings = FALSE)

# Load location combined data
OGlcNAc_location_combined <- read_csv(
  paste0(source_file_path, 'subcellular_location/OGlcNAc_protein_location_combined.csv')
)
cat("Loaded: OGlcNAc_protein_location_combined.csv\n")

# Load normalized data for each cell type
OGlcNAc_protein_norm_HEK293T <- read_csv(
  paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_HEK293T.csv')
)
OGlcNAc_protein_norm_HepG2 <- read_csv(
  paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_HepG2.csv')
)

# Function to extract and calculate log2FC for selected proteins
extract_log2FC_data <- function(norm_data, cell_name, protein_ids) {
  norm_data |>
    filter(Protein.ID %in% protein_ids) |>
    dplyr::select(Protein.ID, Gene,
                  Intensity.Tuni_1_sl_tmm, Intensity.Tuni_2_sl_tmm, Intensity.Tuni_3_sl_tmm,
                  Intensity.Ctrl_4_sl_tmm, Intensity.Ctrl_5_sl_tmm, Intensity.Ctrl_6_sl_tmm) |>
    mutate(
      mean_Ctrl = (Intensity.Ctrl_4_sl_tmm + Intensity.Ctrl_5_sl_tmm + Intensity.Ctrl_6_sl_tmm) / 3,
      log2FC_rep1 = log2(Intensity.Tuni_1_sl_tmm / mean_Ctrl),
      log2FC_rep2 = log2(Intensity.Tuni_2_sl_tmm / mean_Ctrl),
      log2FC_rep3 = log2(Intensity.Tuni_3_sl_tmm / mean_Ctrl),
      cell = cell_name
    ) |>
    dplyr::select(Protein.ID, Gene, cell, log2FC_rep1, log2FC_rep2, log2FC_rep3) |>
    pivot_longer(
      cols = starts_with("log2FC"),
      names_to = "replicate",
      values_to = "log2FC"
    ) |>
    mutate(replicate = str_replace(replicate, "log2FC_rep", ""))
}

# Extract and save location-specific protein lists
OGlcNAc_PM_HEK293T <- OGlcNAc_location_combined |>
  filter(cell == "HEK293T", Location == "Plasma membrane")
write_csv(OGlcNAc_PM_HEK293T, paste0(source_file_path, 'subcellular_location/OGlcNAc_PM_HEK293T.csv'))

OGlcNAc_ER_HepG2 <- OGlcNAc_location_combined |>
  filter(cell == "HepG2", Location == "Endoplasmic reticulum")
write_csv(OGlcNAc_ER_HepG2, paste0(source_file_path, 'subcellular_location/OGlcNAc_ER_HepG2.csv'))

OGlcNAc_Mito_HepG2 <- OGlcNAc_location_combined |>
  filter(cell == "HepG2", Location == "Mitochondria")
write_csv(OGlcNAc_Mito_HepG2, paste0(source_file_path, 'subcellular_location/OGlcNAc_Mito_HepG2.csv'))

cat("Saved location-specific protein lists\n")

# Figure 5C - Plasma Membrane proteins (HEK293T)
selected_proteins_PM <- c("P11166", "O00161", "Q9Y5M8")

log2FC_data_PM <- extract_log2FC_data(OGlcNAc_protein_norm_HEK293T, "HEK293T", selected_proteins_PM) |>
  filter(!is.na(log2FC) & is.finite(log2FC))

cat("\n========== Figure 5C: PM O-GlcNAc Proteins (HEK293T) ==========\n")
print(log2FC_data_PM |> group_by(Gene) |> summarise(n = n(), mean_log2FC = mean(log2FC), .groups = "drop"))

gene_labels_PM <- log2FC_data_PM |> dplyr::select(Protein.ID, Gene) |> distinct()
log2FC_data_PM <- log2FC_data_PM |> mutate(Gene = factor(Gene, levels = gene_labels_PM$Gene))

figure5C <- log2FC_data_PM |>
  ggplot(aes(x = log2FC, y = Gene)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "solid", linewidth = 0.5) +
  geom_jitter(height = 0.15, size = 2, alpha = 0.8, color = colors_cell["HEK293T"]) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.3, color = "black") +
  scale_x_continuous(limits = c(-2, 2)) +
  labs(x = expression(log[2](Tuni/Ctrl)), y = "") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 9),
    axis.text.x = element_text(color = "black", size = 9),
    axis.text.y = element_text(color = "black", size = 9),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = paste0(figure_file_path, 'Figure5/Figure5C.pdf'),
  plot = figure5C,
  height = 1, width = 2, units = 'in'
)
cat("Figure 5C saved\n")

# Figure 5D - ER proteins (HepG2)
selected_proteins_ER <- c("P14314", "Q15084", "Q8NBS9")

log2FC_data_ER <- extract_log2FC_data(OGlcNAc_protein_norm_HepG2, "HepG2", selected_proteins_ER) |>
  filter(!is.na(log2FC) & is.finite(log2FC))

cat("\n========== Figure 5D: ER O-GlcNAc Proteins (HepG2) ==========\n")
print(log2FC_data_ER |> group_by(Gene) |> summarise(n = n(), mean_log2FC = mean(log2FC), .groups = "drop"))

gene_labels_ER <- log2FC_data_ER |> dplyr::select(Protein.ID, Gene) |> distinct()
log2FC_data_ER <- log2FC_data_ER |> mutate(Gene = factor(Gene, levels = gene_labels_ER$Gene))

figure5D <- log2FC_data_ER |>
  ggplot(aes(x = log2FC, y = Gene)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "solid", linewidth = 0.5) +
  geom_jitter(height = 0.15, size = 2, alpha = 0.8, color = colors_cell["HepG2"]) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.3, color = "black") +
  scale_x_continuous(limits = c(-2.15, 2)) +
  labs(x = expression(log[2](Tuni/Ctrl)), y = "") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 9),
    axis.text.x = element_text(color = "black", size = 9),
    axis.text.y = element_text(color = "black", size = 9),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = paste0(figure_file_path, 'Figure5/Figure5D.pdf'),
  plot = figure5D,
  height = 1, width = 2, units = 'in'
)
cat("Figure 5D saved\n")

# Figure 5E - Mitochondria proteins (HepG2)
selected_proteins_Mito <- c("Q9H479", "P28331", "Q92667")

log2FC_data_Mito <- extract_log2FC_data(OGlcNAc_protein_norm_HepG2, "HepG2", selected_proteins_Mito) |>
  filter(!is.na(log2FC) & is.finite(log2FC))

cat("\n========== Figure 5E: Mitochondria O-GlcNAc Proteins (HepG2) ==========\n")
print(log2FC_data_Mito |> group_by(Gene) |> summarise(n = n(), mean_log2FC = mean(log2FC), .groups = "drop"))

gene_labels_Mito <- log2FC_data_Mito |> dplyr::select(Protein.ID, Gene) |> distinct()
log2FC_data_Mito <- log2FC_data_Mito |> mutate(Gene = factor(Gene, levels = gene_labels_Mito$Gene))

figure5E <- log2FC_data_Mito |>
  ggplot(aes(x = log2FC, y = Gene)) +
  geom_vline(xintercept = 0, color = "grey50", linetype = "solid", linewidth = 0.5) +
  geom_jitter(height = 0.15, size = 2, alpha = 0.8, color = colors_cell["HepG2"]) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.3, color = "black") +
  scale_x_continuous(limits = c(-2, 2)) +
  labs(x = expression(log[2](Tuni/Ctrl)), y = "") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 9),
    axis.text.x = element_text(color = "black", size = 9),
    axis.text.y = element_text(color = "black", size = 9),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = paste0(figure_file_path, 'Figure5/Figure5E.pdf'),
  plot = figure5E,
  height = 1, width = 2, units = 'in'
)
cat("Figure 5E saved\n")

cat("\n=== Figure 5C-5E complete ===\n")
