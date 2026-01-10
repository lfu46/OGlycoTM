# Figure 3A Alternative: Combined UMAP Analysis
# Uses Python UMAP via reticulate

library(tidyverse)
library(reticulate)

# Source data paths and colors
source('data_source.R')

# =============================================================================
# SETUP: Configure Python environment (run once)
# =============================================================================
# To create the conda environment, run in terminal:
# conda env create -f umap_env.yml

# Tell reticulate where conda is located
options(reticulate.conda_binary = "/opt/anaconda3/bin/conda")

# Use the umap_env conda environment (using full path)
use_python("/opt/anaconda3/envs/umap_env/bin/python", required = TRUE)

# Source the Python UMAP script
source_python("umap_analysis.py")

# =============================================================================
# Load normalized protein data
# =============================================================================
OGlcNAc_protein_norm_HEK293T <- read_csv(
  paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_HEK293T.csv')
)
OGlcNAc_protein_norm_HepG2 <- read_csv(
  paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_HepG2.csv')
)
OGlcNAc_protein_norm_Jurkat <- read_csv(
  paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_Jurkat.csv')
)

# =============================================================================
# Prepare combined data matrix
# =============================================================================

# Define intensity columns
intensity_cols <- c(
  "Intensity.Tuni_1_sl_tmm", "Intensity.Tuni_2_sl_tmm", "Intensity.Tuni_3_sl_tmm",
  "Intensity.Ctrl_4_sl_tmm", "Intensity.Ctrl_5_sl_tmm", "Intensity.Ctrl_6_sl_tmm"
)

# Function to prepare data for one cell type
prepare_cell_data <- function(data, cell_type) {
  # Assume first column is protein identifier (adjust if needed)
  protein_col <- colnames(data)[1]

  data %>%
    select(all_of(c(protein_col, intensity_cols))) %>%
    rename(Protein = all_of(protein_col)) %>%
    mutate(CellType = cell_type)
}

# Prepare each cell type
data_HEK293T <- prepare_cell_data(OGlcNAc_protein_norm_HEK293T, "HEK293T")
data_HepG2 <- prepare_cell_data(OGlcNAc_protein_norm_HepG2, "HepG2")
data_Jurkat <- prepare_cell_data(OGlcNAc_protein_norm_Jurkat, "Jurkat")

# Find proteins common to all three cell types
common_proteins <- Reduce(intersect, list(
  data_HEK293T$Protein,
  data_HepG2$Protein,
  data_Jurkat$Protein
))

cat("Number of proteins in each cell type:\n")
cat("  HEK293T:", nrow(data_HEK293T), "\n")
cat("  HepG2:", nrow(data_HepG2), "\n")
cat("  Jurkat:", nrow(data_Jurkat), "\n")
cat("Common proteins across all cell types:", length(common_proteins), "\n")

# Filter to common proteins and combine
data_combined <- bind_rows(
  data_HEK293T %>% filter(Protein %in% common_proteins),
  data_HepG2 %>% filter(Protein %in% common_proteins),
  data_Jurkat %>% filter(Protein %in% common_proteins)
)

# =============================================================================
# Reshape data: samples as rows, proteins as columns
# =============================================================================

# Create sample metadata
sample_info <- tibble(
  Sample = c(
    paste0("HEK293T_Tuni_", 1:3), paste0("HEK293T_Ctrl_", 1:3),
    paste0("HepG2_Tuni_", 1:3), paste0("HepG2_Ctrl_", 1:3),
    paste0("Jurkat_Tuni_", 1:3), paste0("Jurkat_Ctrl_", 1:3)
  ),
  CellType = rep(c("HEK293T", "HepG2", "Jurkat"), each = 6),
  Condition = rep(rep(c("Tuni", "Ctrl"), each = 3), 3),
  OriginalCol = rep(intensity_cols, 3)
)

# Create wide matrix: rows = samples, columns = proteins
create_sample_matrix <- function(cell_data, cell_type) {
  cell_data %>%
    filter(CellType == cell_type) %>%
    select(Protein, all_of(intensity_cols)) %>%
    pivot_longer(
      cols = all_of(intensity_cols),
      names_to = "SampleCol",
      values_to = "Intensity"
    ) %>%
    mutate(
      Condition = ifelse(str_detect(SampleCol, "Tuni"), "Tuni", "Ctrl"),
      Replicate = str_extract(SampleCol, "\\d"),
      Sample = paste0(cell_type, "_", Condition, "_", Replicate)
    ) %>%
    select(Sample, Protein, Intensity) %>%
    pivot_wider(names_from = Protein, values_from = Intensity)
}

# Create matrices for each cell type
matrix_HEK293T <- create_sample_matrix(data_combined, "HEK293T")
matrix_HepG2 <- create_sample_matrix(data_combined, "HepG2")
matrix_Jurkat <- create_sample_matrix(data_combined, "Jurkat")

# Combine all samples
intensity_matrix <- bind_rows(matrix_HEK293T, matrix_HepG2, matrix_Jurkat)

# Extract sample names and prepare for UMAP
sample_names <- intensity_matrix$Sample
intensity_data <- intensity_matrix %>% select(-Sample)

# Handle NA/zero values: log2 transform and remove problematic proteins
intensity_data <- log2(intensity_data + 1)

# Remove proteins with any NA values
intensity_data <- intensity_data %>%
  select(where(~ !any(is.na(.))))

cat("Final matrix dimensions:", nrow(intensity_data), "samples x",
    ncol(intensity_data), "proteins\n")

# =============================================================================
# Run UMAP via Python
# =============================================================================

# Prepare metadata vectors
cell_types <- str_extract(sample_names, "^[^_]+")
conditions <- str_extract(sample_names, "(?<=_)[^_]+(?=_)")

# Run UMAP (adjust n_neighbors for small sample size)
umap_result <- run_umap_with_metadata(
  intensity_df = intensity_data,
  sample_names = sample_names,
  cell_types = cell_types,
  conditions = conditions,
  n_neighbors = 5L,  # Use integer for Python
  min_dist = 0.3,
  random_state = 42L
)

# Convert to tibble for ggplot
umap_df <- as_tibble(umap_result)

# =============================================================================
# Create UMAP visualization
# =============================================================================

# Define colors and shapes
colors_treatment <- c("Ctrl" = "#8491B4", "Tuni" = "#E64B35")

# Plot with cell type as color, condition as shape
p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = CellType, shape = Condition), size = 4, stroke = 1.2) +
  scale_color_manual(values = colors_cell) +
  scale_shape_manual(values = c("Ctrl" = 16, "Tuni" = 17)) +  # circle, triangle
  labs(
    x = "UMAP1",
    y = "UMAP2",
    title = "O-GlcNAc Protein UMAP: Tuni vs Ctrl across Cell Types",
    color = "Cell Type",
    shape = "Condition"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

print(p_umap)

# Alternative: Color by condition, shape by cell type
p_umap_alt <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = Condition, shape = CellType), size = 4, stroke = 1.2) +
  scale_color_manual(values = colors_treatment) +
  scale_shape_manual(values = c("HEK293T" = 16, "HepG2" = 17, "Jurkat" = 15)) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    title = "O-GlcNAc Protein UMAP: Tuni vs Ctrl across Cell Types",
    color = "Condition",
    shape = "Cell Type"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

print(p_umap_alt)

# =============================================================================
# Save plots
# =============================================================================

ggsave(
  paste0(figure_file_path, "Figure3A_UMAP_combined.pdf"),
  p_umap,
  width = 7, height = 5
)

ggsave(
  paste0(figure_file_path, "Figure3A_UMAP_combined_alt.pdf"),
  p_umap_alt,
  width = 7, height = 5
)

cat("\nUMAP plots saved to:", figure_file_path, "\n")

# =============================================================================
# Optional: Print UMAP coordinates for inspection
# =============================================================================
cat("\nUMAP Coordinates:\n")
print(umap_df)
