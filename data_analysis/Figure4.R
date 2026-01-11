# Figure 4: Circular Heatmap Visualization of O-GlcNAc Protein Differential Expression

library(tidyverse)
library(circlize)
library(scales)
library(ComplexHeatmap)
library(gridBase)

# Source data paths, colors, and differential analysis results
source('data_source_DE.R')

# =============================================================================
# Figure 4A - Circular Heatmap with Intralinks
# =============================================================================
# This figure shows all O-GlcNAc proteins across three cell types (HEK293T, HepG2, Jurkat)
# with multiple layers showing:
# - Category (up/median/down regulation)
# - Z-score normalized TMT abundance (6 channels)
# - Cell type
# - Adjusted p-value
# Intralinks connect proteins that are significantly down in one cell type
# but significantly up in another cell type

# Load normalized O-GlcNAc protein data
OGlcNAc_protein_norm_HEK293T <- read_csv(
  paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_HEK293T.csv')
)
OGlcNAc_protein_norm_HepG2 <- read_csv(
  paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_HepG2.csv')
)
OGlcNAc_protein_norm_Jurkat <- read_csv(
  paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_Jurkat.csv')
)

# Function to categorize and prepare data for one cell type
# Categories: up (logFC > 0.5, adj.P.Val < 0.05), down (logFC < -0.5, adj.P.Val < 0.05), median (rest)
prepare_cell_data <- function(de_data, norm_data, cell_name) {
  up_proteins <- de_data |>
    filter(logFC > 0.5, adj.P.Val < 0.05) |>
    dplyr::select(Protein.ID, logFC, adj.P.Val) |>
    mutate(category = "up", cell = cell_name) |>
    arrange(adj.P.Val) |>
    left_join(norm_data |> dplyr::select(Protein.ID, ends_with("_sl_tmm")), by = "Protein.ID")

  down_proteins <- de_data |>
    filter(logFC < -0.5, adj.P.Val < 0.05) |>
    dplyr::select(Protein.ID, logFC, adj.P.Val) |>
    mutate(category = "down", cell = cell_name) |>
    arrange(adj.P.Val) |>
    left_join(norm_data |> dplyr::select(Protein.ID, ends_with("_sl_tmm")), by = "Protein.ID")

  median_proteins <- de_data |>
    filter(!(Protein.ID %in% up_proteins$Protein.ID)) |>
    filter(!(Protein.ID %in% down_proteins$Protein.ID)) |>
    dplyr::select(Protein.ID, logFC, adj.P.Val) |>
    mutate(category = "median", cell = cell_name) |>
    arrange(adj.P.Val) |>
    left_join(norm_data |> dplyr::select(Protein.ID, ends_with("_sl_tmm")), by = "Protein.ID")

  list(up = up_proteins, median = median_proteins, down = down_proteins)
}

# Prepare data for each cell type
HepG2_data <- prepare_cell_data(OGlcNAc_protein_DE_HepG2, OGlcNAc_protein_norm_HepG2, "HepG2")
HEK293T_data <- prepare_cell_data(OGlcNAc_protein_DE_HEK293T, OGlcNAc_protein_norm_HEK293T, "HEK293T")
Jurkat_data <- prepare_cell_data(OGlcNAc_protein_DE_Jurkat, OGlcNAc_protein_norm_Jurkat, "Jurkat")

# Print counts
cat("HepG2: up =", nrow(HepG2_data$up), ", median =", nrow(HepG2_data$median), ", down =", nrow(HepG2_data$down), "\n")
cat("HEK293T: up =", nrow(HEK293T_data$up), ", median =", nrow(HEK293T_data$median), ", down =", nrow(HEK293T_data$down), "\n")
cat("Jurkat: up =", nrow(Jurkat_data$up), ", median =", nrow(Jurkat_data$median), ", down =", nrow(Jurkat_data$down), "\n")

# Combine data: HepG2 (up, median, down) -> HEK293T (up, median, down) -> Jurkat (up, median, down)
OGlcNAc_protein_combined <- bind_rows(
  HepG2_data$up, HepG2_data$median, HepG2_data$down,
  HEK293T_data$up, HEK293T_data$median, HEK293T_data$down,
  Jurkat_data$up, Jurkat_data$median, Jurkat_data$down
)

# Row-wise min-max scaling to [-1, 1] for Z-score visualization
intensity_cols <- names(OGlcNAc_protein_combined)[str_detect(names(OGlcNAc_protein_combined), "_sl_tmm$")]
OGlcNAc_protein_combined <- OGlcNAc_protein_combined |>
  rowwise() |>
  mutate(
    max_value = max(c_across(all_of(intensity_cols)), na.rm = TRUE),
    min_value = min(c_across(all_of(intensity_cols)), na.rm = TRUE)
  ) |>
  ungroup()

for (col in intensity_cols) {
  scaled_col <- paste0("scaled_", gsub("Intensity\\.", "", col))
  OGlcNAc_protein_combined <- OGlcNAc_protein_combined |>
    mutate(!!scaled_col := (.data[[col]] - min_value) * 2 / (max_value - min_value) - 1)
}

# Add cell_line index (1 = HEK293T, 2 = HepG2, 3 = Jurkat)
OGlcNAc_protein_combined <- OGlcNAc_protein_combined |>
  mutate(cell_line = case_when(cell == "HEK293T" ~ 1, cell == "HepG2" ~ 2, cell == "Jurkat" ~ 3))

# Save combined data
dir.create(paste0(source_file_path, 'circular_heatmap'), showWarnings = FALSE)
write_csv(OGlcNAc_protein_combined, paste0(source_file_path, 'circular_heatmap/OGlcNAc_protein_circular_heatmap_data.csv'))

# Prepare data matrices for circular heatmap
scaled_cols <- names(OGlcNAc_protein_combined)[str_detect(names(OGlcNAc_protein_combined), "^scaled_")]
mat_sl_tmm <- data.matrix(OGlcNAc_protein_combined |> dplyr::select(all_of(scaled_cols)))
cell <- OGlcNAc_protein_combined$cell
cell_line <- OGlcNAc_protein_combined$cell_line
adjpVal <- OGlcNAc_protein_combined$adj.P.Val
category <- OGlcNAc_protein_combined$category

# Function to find overlapping proteins and their indices for intralinks
find_link_indices <- function(down_data, up_data) {
  overlap <- semi_join(down_data, up_data, by = "Protein.ID")
  if (nrow(overlap) == 0) return(tibble(down_index = integer(), up_index = integer()))
  down_indices <- which(down_data$Protein.ID %in% overlap$Protein.ID)
  up_indices <- which(up_data$Protein.ID %in% overlap$Protein.ID)
  down_ids <- down_data$Protein.ID[down_indices]
  up_ids <- up_data$Protein.ID[up_indices]
  matched <- tibble(Protein.ID = down_ids) |>
    mutate(down_index = down_indices) |>
    left_join(tibble(Protein.ID = up_ids, up_index = up_indices), by = "Protein.ID") |>
    filter(!is.na(up_index))
  matched |> dplyr::select(down_index, up_index)
}

# Calculate offsets (position of "down" section within each cell sector)
offset_HepG2 <- nrow(HepG2_data$up) + nrow(HepG2_data$median)
offset_HEK293T <- nrow(HEK293T_data$up) + nrow(HEK293T_data$median)
offset_Jurkat <- nrow(Jurkat_data$up) + nrow(Jurkat_data$median)
cat("\nOffsets: HepG2 =", offset_HepG2, ", HEK293T =", offset_HEK293T, ", Jurkat =", offset_Jurkat, "\n")

# Find all link pairs
link_downHepG2_upHEK293T <- find_link_indices(HepG2_data$down, HEK293T_data$up)
link_downHepG2_upJurkat <- find_link_indices(HepG2_data$down, Jurkat_data$up)
link_downHEK293T_upHepG2 <- find_link_indices(HEK293T_data$down, HepG2_data$up)
link_downHEK293T_upJurkat <- find_link_indices(HEK293T_data$down, Jurkat_data$up)
link_downJurkat_upHEK293T <- find_link_indices(Jurkat_data$down, HEK293T_data$up)
link_downJurkat_upHepG2 <- find_link_indices(Jurkat_data$down, HepG2_data$up)

cat("\nLink counts:\n")
cat("downHepG2 -> upHEK293T:", nrow(link_downHepG2_upHEK293T), "\n")
cat("downHepG2 -> upJurkat:", nrow(link_downHepG2_upJurkat), "\n")
cat("downHEK293T -> upHepG2:", nrow(link_downHEK293T_upHepG2), "\n")
cat("downHEK293T -> upJurkat:", nrow(link_downHEK293T_upJurkat), "\n")
cat("downJurkat -> upHEK293T:", nrow(link_downJurkat_upHEK293T), "\n")
cat("downJurkat -> upHepG2:", nrow(link_downJurkat_upHepG2), "\n")

# Define colors
col_category <- c("up" = "#FFCC00", "median" = "gray80", "down" = "#7B68EE")
col_mat <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))
col_cell <- colors_cell
col_adjpval <- colorRamp2(c(0, 0.001, 0.01, 0.05, 1), c("#006400", "#228B22", "#32CD32", "#90EE90", "white"))

# Category counts and positions for labels
n_up_HepG2 <- nrow(HepG2_data$up); n_median_HepG2 <- nrow(HepG2_data$median); n_down_HepG2 <- nrow(HepG2_data$down)
n_up_HEK293T <- nrow(HEK293T_data$up); n_median_HEK293T <- nrow(HEK293T_data$median); n_down_HEK293T <- nrow(HEK293T_data$down)
n_up_Jurkat <- nrow(Jurkat_data$up); n_median_Jurkat <- nrow(Jurkat_data$median); n_down_Jurkat <- nrow(Jurkat_data$down)

pos_up_HepG2 <- n_up_HepG2 / 2; pos_median_HepG2 <- n_up_HepG2 + n_median_HepG2 / 2; pos_down_HepG2 <- n_up_HepG2 + n_median_HepG2 + n_down_HepG2 / 2
pos_up_HEK293T <- n_up_HEK293T / 2; pos_median_HEK293T <- n_up_HEK293T + n_median_HEK293T / 2; pos_down_HEK293T <- n_up_HEK293T + n_median_HEK293T + n_down_HEK293T / 2
pos_up_Jurkat <- n_up_Jurkat / 2; pos_median_Jurkat <- n_up_Jurkat + n_median_Jurkat / 2; pos_down_Jurkat <- n_up_Jurkat + n_median_Jurkat + n_down_Jurkat / 2

# Circular heatmap function
circlize_plot <- function() {
  # Track 1: Category (outermost)
  circos.heatmap(category, split = cell_line, col = col_category, track.height = 0.03)

  # Add category count labels outside the heatmap
  # Note: y > 1 places text outside track bounds, which generates harmless warnings
  circos.text(pos_up_HEK293T, 2.5, n_up_HEK293T, sector.index = "1", col = col_cell["HEK293T"], cex = 0.6, font = 2)
  circos.text(pos_median_HEK293T, 2.5, n_median_HEK293T, sector.index = "1", col = "black", cex = 0.6, font = 2)
  circos.text(pos_down_HEK293T, 2.5, n_down_HEK293T, sector.index = "1", col = col_cell["HEK293T"], cex = 0.6, font = 2)

  circos.text(pos_up_HepG2, 2.5, n_up_HepG2, sector.index = "2", col = col_cell["HepG2"], cex = 0.6, font = 2)
  circos.text(pos_median_HepG2, 2.5, n_median_HepG2, sector.index = "2", col = "black", cex = 0.6, font = 2)
  circos.text(pos_down_HepG2, 2.5, n_down_HepG2, sector.index = "2", col = col_cell["HepG2"], cex = 0.6, font = 2)

  circos.text(pos_up_Jurkat, 2.5, n_up_Jurkat, sector.index = "3", col = col_cell["Jurkat"], cex = 0.6, font = 2)
  circos.text(pos_median_Jurkat, 2.5, n_median_Jurkat, sector.index = "3", col = "black", cex = 0.6, font = 2)
  circos.text(pos_down_Jurkat, 2.5, n_down_Jurkat, sector.index = "3", col = col_cell["Jurkat"], cex = 0.6, font = 2)

  # Track 2: Scaled TMT abundance
  circos.heatmap(mat_sl_tmm, col = col_mat, track.height = 0.12)

  # Track 3: Cell type
  circos.heatmap(cell, col = col_cell, track.height = 0.03)

  # Track 4: Adjusted p-value (innermost)
  circos.heatmap(adjpVal, col = col_adjpval, track.height = 0.03)

  # Add intralinks (colored by source cell type where protein is down)
  if (nrow(link_downHepG2_upHEK293T) > 0) {
    for (i in seq_len(nrow(link_downHepG2_upHEK293T))) {
      circos.link(2, link_downHepG2_upHEK293T$down_index[i] + offset_HepG2 - 0.5,
                  1, link_downHepG2_upHEK293T$up_index[i] - 0.5, col = alpha(col_cell["HepG2"], 0.6), lwd = 2)
    }
  }
  if (nrow(link_downHepG2_upJurkat) > 0) {
    for (i in seq_len(nrow(link_downHepG2_upJurkat))) {
      circos.link(2, link_downHepG2_upJurkat$down_index[i] + offset_HepG2 - 0.5,
                  3, link_downHepG2_upJurkat$up_index[i] - 0.5, col = alpha(col_cell["HepG2"], 0.6), lwd = 2)
    }
  }
  if (nrow(link_downHEK293T_upHepG2) > 0) {
    for (i in seq_len(nrow(link_downHEK293T_upHepG2))) {
      circos.link(1, link_downHEK293T_upHepG2$down_index[i] + offset_HEK293T - 0.5,
                  2, link_downHEK293T_upHepG2$up_index[i] - 0.5, col = alpha(col_cell["HEK293T"], 0.6), lwd = 2)
    }
  }
  if (nrow(link_downHEK293T_upJurkat) > 0) {
    for (i in seq_len(nrow(link_downHEK293T_upJurkat))) {
      circos.link(1, link_downHEK293T_upJurkat$down_index[i] + offset_HEK293T - 0.5,
                  3, link_downHEK293T_upJurkat$up_index[i] - 0.5, col = alpha(col_cell["HEK293T"], 0.6), lwd = 2)
    }
  }
  if (nrow(link_downJurkat_upHEK293T) > 0) {
    for (i in seq_len(nrow(link_downJurkat_upHEK293T))) {
      circos.link(3, link_downJurkat_upHEK293T$down_index[i] + offset_Jurkat - 0.5,
                  1, link_downJurkat_upHEK293T$up_index[i] - 0.5, col = alpha(col_cell["Jurkat"], 0.6), lwd = 2)
    }
  }
  if (nrow(link_downJurkat_upHepG2) > 0) {
    for (i in seq_len(nrow(link_downJurkat_upHepG2))) {
      circos.link(3, link_downJurkat_upHepG2$down_index[i] + offset_Jurkat - 0.5,
                  2, link_downJurkat_upHepG2$up_index[i] - 0.5, col = alpha(col_cell["Jurkat"], 0.6), lwd = 2)
    }
  }

  circos.clear()
}

# Generate legends (smaller size)
lgd_mat <- Legend(title = "Z-score", col_fun = col_mat, title_gp = gpar(fontsize = 6, fontface = "bold"),
                  labels_gp = gpar(fontsize = 5), grid_height = unit(3, "mm"), grid_width = unit(3, "mm"), legend_height = unit(12, "mm"))
lgd_cell <- Legend(title = "Cell", at = names(col_cell), legend_gp = gpar(fill = col_cell),
                   title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 5),
                   grid_height = unit(3, "mm"), grid_width = unit(3, "mm"))
lgd_category <- Legend(title = "Category", at = names(col_category), legend_gp = gpar(fill = col_category),
                       title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 5),
                       grid_height = unit(3, "mm"), grid_width = unit(3, "mm"))
lgd_adjpval <- Legend(title = "adj.P.Value", col_fun = col_adjpval, at = c(0, 0.001, 0.01, 0.05, 1),
                      title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 5),
                      grid_height = unit(3, "mm"), grid_width = unit(3, "mm"), legend_height = unit(12, "mm"))

# Create and save the circular heatmap
dir.create(paste0(figure_file_path, "Figure4"), showWarnings = FALSE)

pdf(file = paste0(figure_file_path, "Figure4/Figure4A.pdf"), width = 5, height = 4)

plot.new()
circle_size <- unit(0.95, "snpc")
pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size, just = c("left", "center")))
par(omi = gridOMI(), new = TRUE)
circlize_plot()
upViewport()

h <- dev.size()[2]
lgd_list <- packLegend(lgd_mat, lgd_cell, lgd_category, lgd_adjpval, max_height = unit(0.98 * h, "inch"), gap = unit(1.5, "mm"))
draw(lgd_list, x = unit(0.78, "npc"), just = "left")

dev.off()

cat("\nFigure 4A saved to:", figure_file_path, "Figure4/\n")

# =============================================================================
# Figure 4B - Oppositely Regulated Proteins List
# =============================================================================
# Extract proteins that are significantly upregulated in one cell type
# but significantly downregulated in another cell type

# Get up and down regulated proteins for each cell type
up_HEK293T <- OGlcNAc_protein_DE_HEK293T |> filter(logFC > 0.5, adj.P.Val < 0.05) |> dplyr::select(Protein.ID)
up_HepG2 <- OGlcNAc_protein_DE_HepG2 |> filter(logFC > 0.5, adj.P.Val < 0.05) |> dplyr::select(Protein.ID)
up_Jurkat <- OGlcNAc_protein_DE_Jurkat |> filter(logFC > 0.5, adj.P.Val < 0.05) |> dplyr::select(Protein.ID)

down_HEK293T <- OGlcNAc_protein_DE_HEK293T |> filter(logFC < -0.5, adj.P.Val < 0.05) |> dplyr::select(Protein.ID)
down_HepG2 <- OGlcNAc_protein_DE_HepG2 |> filter(logFC < -0.5, adj.P.Val < 0.05) |> dplyr::select(Protein.ID)
down_Jurkat <- OGlcNAc_protein_DE_Jurkat |> filter(logFC < -0.5, adj.P.Val < 0.05) |> dplyr::select(Protein.ID)

# Find all oppositely regulated proteins (union of all combinations)
oppositely_regulated <- bind_rows(
  # Down in HepG2, Up in HEK293T
  inner_join(down_HepG2, up_HEK293T, by = "Protein.ID") |> mutate(pattern = "down_HepG2_up_HEK293T"),
  # Down in HepG2, Up in Jurkat
  inner_join(down_HepG2, up_Jurkat, by = "Protein.ID") |> mutate(pattern = "down_HepG2_up_Jurkat"),
  # Down in HEK293T, Up in HepG2
  inner_join(down_HEK293T, up_HepG2, by = "Protein.ID") |> mutate(pattern = "down_HEK293T_up_HepG2"),
  # Down in HEK293T, Up in Jurkat
  inner_join(down_HEK293T, up_Jurkat, by = "Protein.ID") |> mutate(pattern = "down_HEK293T_up_Jurkat"),
  # Down in Jurkat, Up in HEK293T
  inner_join(down_Jurkat, up_HEK293T, by = "Protein.ID") |> mutate(pattern = "down_Jurkat_up_HEK293T"),
  # Down in Jurkat, Up in HepG2
  inner_join(down_Jurkat, up_HepG2, by = "Protein.ID") |> mutate(pattern = "down_Jurkat_up_HepG2")
)

# Get unique proteins
unique_proteins <- unique(oppositely_regulated$Protein.ID)
cat("\nTotal unique oppositely regulated proteins:", length(unique_proteins), "\n")

# Create comprehensive table with fold changes from all cell types
oppositely_regulated_table <- tibble(Protein.ID = unique_proteins) |>
  left_join(
    OGlcNAc_protein_DE_HEK293T |> dplyr::select(Protein.ID, logFC_HEK293T = logFC, adjPVal_HEK293T = adj.P.Val),
    by = "Protein.ID"
  ) |>
  left_join(
    OGlcNAc_protein_DE_HepG2 |> dplyr::select(Protein.ID, logFC_HepG2 = logFC, adjPVal_HepG2 = adj.P.Val),
    by = "Protein.ID"
  ) |>
  left_join(
    OGlcNAc_protein_DE_Jurkat |> dplyr::select(Protein.ID, logFC_Jurkat = logFC, adjPVal_Jurkat = adj.P.Val),
    by = "Protein.ID"
  )

# Add regulation status for each cell type
oppositely_regulated_table <- oppositely_regulated_table |>
  mutate(
    status_HEK293T = case_when(
      logFC_HEK293T > 0.5 & adjPVal_HEK293T < 0.05 ~ "Up",
      logFC_HEK293T < -0.5 & adjPVal_HEK293T < 0.05 ~ "Down",
      TRUE ~ "NS"
    ),
    status_HepG2 = case_when(
      logFC_HepG2 > 0.5 & adjPVal_HepG2 < 0.05 ~ "Up",
      logFC_HepG2 < -0.5 & adjPVal_HepG2 < 0.05 ~ "Down",
      TRUE ~ "NS"
    ),
    status_Jurkat = case_when(
      logFC_Jurkat > 0.5 & adjPVal_Jurkat < 0.05 ~ "Up",
      logFC_Jurkat < -0.5 & adjPVal_Jurkat < 0.05 ~ "Down",
      TRUE ~ "NS"
    )
  ) |>
  arrange(Protein.ID)

# Print summary
cat("\nOppositely regulated proteins by pattern:\n")
oppositely_regulated |>
  group_by(pattern) |>
  summarise(n = n(), .groups = "drop") |>
  print()

# Print the complete table
cat("\nComplete list of oppositely regulated proteins:\n")
print(oppositely_regulated_table, n = Inf)

# Save to CSV
write_csv(oppositely_regulated_table, paste0(source_file_path, 'circular_heatmap/oppositely_regulated_proteins.csv'))
cat("\nOppositely regulated proteins table saved to:", source_file_path, "circular_heatmap/oppositely_regulated_proteins.csv\n")

# =============================================================================
# Figure 4B - Dot Plot for Selected Oppositely Regulated Proteins
# =============================================================================
# Dot plot showing log2(Tuni/Ctrl) for selected proteins across cell types
# Proteins: Q8IX12, P27540, Q9UHY1

# Define proteins of interest
proteins_of_interest <- c("Q8IX12", "P27540", "Q9UHY1")

# Function to extract and calculate fold changes for a cell type
calculate_fc <- function(norm_data, cell_name, proteins) {
  norm_data |>
    filter(Protein.ID %in% proteins) |>
    dplyr::select(Protein.ID, Gene,
                  Intensity.Tuni_1_sl_tmm, Intensity.Tuni_2_sl_tmm, Intensity.Tuni_3_sl_tmm,
                  Intensity.Ctrl_4_sl_tmm, Intensity.Ctrl_5_sl_tmm, Intensity.Ctrl_6_sl_tmm) |>
    rowwise() |>
    mutate(
      mean_Ctrl = mean(c(Intensity.Ctrl_4_sl_tmm, Intensity.Ctrl_5_sl_tmm, Intensity.Ctrl_6_sl_tmm), na.rm = TRUE),
      log2FC_rep1 = log2(Intensity.Tuni_1_sl_tmm / mean_Ctrl),
      log2FC_rep2 = log2(Intensity.Tuni_2_sl_tmm / mean_Ctrl),
      log2FC_rep3 = log2(Intensity.Tuni_3_sl_tmm / mean_Ctrl)
    ) |>
    ungroup() |>
    dplyr::select(Protein.ID, Gene, log2FC_rep1, log2FC_rep2, log2FC_rep3) |>
    pivot_longer(cols = starts_with("log2FC"), names_to = "replicate", values_to = "log2FC") |>
    mutate(cell = cell_name)
}

# Calculate fold changes for each cell type
fc_HEK293T <- calculate_fc(OGlcNAc_protein_norm_HEK293T, "HEK293T", proteins_of_interest)
fc_HepG2 <- calculate_fc(OGlcNAc_protein_norm_HepG2, "HepG2", proteins_of_interest)
fc_Jurkat <- calculate_fc(OGlcNAc_protein_norm_Jurkat, "Jurkat", proteins_of_interest)

# Combine all data
fc_combined <- bind_rows(fc_HEK293T, fc_HepG2, fc_Jurkat) |>
  mutate(
    cell = factor(cell, levels = c("HEK293T", "HepG2", "Jurkat")),
    Gene = factor(Gene)
  )

# Check data
cat("\nFold change data for selected proteins:\n")
print(fc_combined)

# Create dot plot
figure4B <- fc_combined |>
  ggplot(aes(x = cell, y = log2FC, color = cell)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0.5, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -0.5, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_point(size = 2, position = position_jitter(width = 0.1, seed = 42)) +
  scale_color_manual(values = colors_cell) +
  facet_wrap(~ Gene, nrow = 1) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 7),
    axis.text.x = element_text(color = "black", size = 7, angle = 30, hjust = 1),
    axis.text.y = element_text(color = "black", size = 7),
    strip.text = element_text(size = 7, face = "bold"),
    strip.background = element_blank(),
    legend.position = "none"
  )

ggsave(
  filename = paste0(figure_file_path, 'Figure4/Figure4B.pdf'),
  plot = figure4B,
  height = 1.8, width = 3.5, units = 'in'
)

cat("\nFigure 4B saved to:", figure_file_path, "Figure4/Figure4B.pdf\n")

# =============================================================================
# Figure 4C - GO Enrichment Analysis for HEK293T
# =============================================================================
# Functional enrichment analysis for up and down regulated O-GlcNAc proteins in HEK293T

library(clusterProfiler)
library(org.Hs.eg.db)

# Create output directory for enrichment results
dir.create(paste0(source_file_path, 'enrichment'), showWarnings = FALSE)

# Define universe: all O-GlcNAc proteins from all three cell types
OGlcNAc_protein_universe <- bind_rows(
  OGlcNAc_protein_DE_HEK293T |> dplyr::select(Protein.ID),
  OGlcNAc_protein_DE_HepG2 |> dplyr::select(Protein.ID),
  OGlcNAc_protein_DE_Jurkat |> dplyr::select(Protein.ID)
) |>
  distinct() |>
  pull(Protein.ID)

cat("\nUniverse size (total O-GlcNAc proteins):", length(OGlcNAc_protein_universe), "\n")

# Upregulated O-GlcNAc proteins in HEK293T
OGlcNAc_up_HEK293T <- OGlcNAc_protein_DE_HEK293T |>
  filter(logFC > 0.5, adj.P.Val < 0.05) |>
  pull(Protein.ID)

# Downregulated O-GlcNAc proteins in HEK293T
OGlcNAc_down_HEK293T <- OGlcNAc_protein_DE_HEK293T |>
  filter(logFC < -0.5, adj.P.Val < 0.05) |>
  pull(Protein.ID)

cat("HEK293T - Up:", length(OGlcNAc_up_HEK293T), ", Down:", length(OGlcNAc_down_HEK293T), "\n")

# GO enrichment for upregulated proteins
OGlcNAc_up_HEK293T_GO <- enrichGO(
  gene = OGlcNAc_up_HEK293T,
  OrgDb = org.Hs.eg.db,
  universe = OGlcNAc_protein_universe,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  OGlcNAc_up_HEK293T_GO@result,
  file = paste0(source_file_path, 'enrichment/OGlcNAc_up_HEK293T_GO.csv')
)

# GO enrichment for downregulated proteins
OGlcNAc_down_HEK293T_GO <- enrichGO(
  gene = OGlcNAc_down_HEK293T,
  OrgDb = org.Hs.eg.db,
  universe = OGlcNAc_protein_universe,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  OGlcNAc_down_HEK293T_GO@result,
  file = paste0(source_file_path, 'enrichment/OGlcNAc_down_HEK293T_GO.csv')
)

cat("HEK293T GO enrichment results saved\n")

# Import saved results for barplot
OGlcNAc_up_HEK293T_GO_result <- read_csv(
  paste0(source_file_path, 'enrichment/OGlcNAc_up_HEK293T_GO.csv')
)
OGlcNAc_down_HEK293T_GO_result <- read_csv(
  paste0(source_file_path, 'enrichment/OGlcNAc_down_HEK293T_GO.csv')
)

# Preview top significant terms
cat("\nHEK293T Up - Top 10 significant terms:\n")
OGlcNAc_up_HEK293T_GO_result |>
  filter(pvalue < 0.05) |>
  arrange(pvalue) |>
  dplyr::select(ONTOLOGY, Description, pvalue, Count) |>
  head(10) |>
  print()

cat("\nHEK293T Down - Top 10 significant terms:\n")
OGlcNAc_down_HEK293T_GO_result |>
  filter(pvalue < 0.05) |>
  arrange(pvalue) |>
  dplyr::select(ONTOLOGY, Description, pvalue, Count) |>
  head(10) |>
  print()

# Barplot for HEK293T (update Description filter after reviewing results)
# TODO: Update the Description list after manual review of enrichment results
figure4C <- OGlcNAc_up_HEK293T_GO_result |>
  filter(
    Description %in% c(
      'ribonucleotide metabolic process', 
      'ribose phosphate metabolic process', 
      'protein folding chaperone', 
      'cytosol', 
      'nucleotide binding'
    )
  ) |>
  ggplot() +
  geom_bar(
    aes(x = -log10(pvalue), y = fct_reorder(Description, -log10(pvalue))),
    stat = "identity", fill = unname(colors_cell["HEK293T"]), width = 0.7
  ) +
  labs(x = expression(-log[10]*"("*paste(italic(P), " Value")*")"), y = "") +
  geom_text(aes(label = Description, x = 0, y = Description), hjust = 0, size = 2.5) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(color = "black", size = 7),
    axis.text.y = element_blank(),
    axis.ticks.length.y = unit(0, "in")
  )

ggsave(
  filename = paste0(figure_file_path, 'Figure4/Figure4C.pdf'),
  plot = figure4C,
  height = 1.5, width = 2, units = 'in'
)

# =============================================================================
# Figure 4D - GO Enrichment Analysis for HepG2
# =============================================================================
# Functional enrichment analysis for up and down regulated O-GlcNAc proteins in HepG2

# Upregulated O-GlcNAc proteins in HepG2
OGlcNAc_up_HepG2 <- OGlcNAc_protein_DE_HepG2 |>
  filter(logFC > 0.5, adj.P.Val < 0.05) |>
  pull(Protein.ID)

# Downregulated O-GlcNAc proteins in HepG2
OGlcNAc_down_HepG2 <- OGlcNAc_protein_DE_HepG2 |>
  filter(logFC < -0.5, adj.P.Val < 0.05) |>
  pull(Protein.ID)

cat("\nHepG2 - Up:", length(OGlcNAc_up_HepG2), ", Down:", length(OGlcNAc_down_HepG2), "\n")

# GO enrichment for upregulated proteins
OGlcNAc_up_HepG2_GO <- enrichGO(
  gene = OGlcNAc_up_HepG2,
  OrgDb = org.Hs.eg.db,
  universe = OGlcNAc_protein_universe,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  OGlcNAc_up_HepG2_GO@result,
  file = paste0(source_file_path, 'enrichment/OGlcNAc_up_HepG2_GO.csv')
)

# GO enrichment for downregulated proteins
OGlcNAc_down_HepG2_GO <- enrichGO(
  gene = OGlcNAc_down_HepG2,
  OrgDb = org.Hs.eg.db,
  universe = OGlcNAc_protein_universe,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  OGlcNAc_down_HepG2_GO@result,
  file = paste0(source_file_path, 'enrichment/OGlcNAc_down_HepG2_GO.csv')
)

cat("HepG2 GO enrichment results saved\n")

# Import saved results for barplot
OGlcNAc_up_HepG2_GO_result <- read_csv(
  paste0(source_file_path, 'enrichment/OGlcNAc_up_HepG2_GO.csv')
)
OGlcNAc_down_HepG2_GO_result <- read_csv(
  paste0(source_file_path, 'enrichment/OGlcNAc_down_HepG2_GO.csv')
)

# Preview top significant terms
cat("\nHepG2 Up - Top 10 significant terms:\n")
OGlcNAc_up_HepG2_GO_result |>
  filter(pvalue < 0.05) |>
  arrange(pvalue) |>
  dplyr::select(ONTOLOGY, Description, pvalue, Count) |>
  head(10) |>
  print()

cat("\nHepG2 Down - Top 10 significant terms:\n")
OGlcNAc_down_HepG2_GO_result |>
  filter(pvalue < 0.05) |>
  arrange(pvalue) |>
  dplyr::select(ONTOLOGY, Description, pvalue, Count) |>
  head(10) |>
  print()

# Barplot for HepG2 (update Description filter after reviewing results)
# TODO: Update the Description list after manual review of enrichment results
figure4D <- OGlcNAc_down_HepG2_GO_result |>
  filter(
    Description %in% c(
      'endoplasmic reticulum lumen',
      'structural constituent of nuclear pore',
      'endoplasmic reticulum',
      'import into nucleus',
      'RNA export from nucleus'
    )
  ) |>
  mutate(Description = str_replace(Description, "endoplasmic reticulum", "ER")) |>
  ggplot() +
  geom_bar(
    aes(x = -log10(pvalue), y = fct_reorder(Description, -log10(pvalue))),
    stat = "identity", fill = unname(colors_cell["HepG2"]), width = 0.7
  ) +
  labs(x = expression(-log[10]*"("*paste(italic(P), " Value")*")"), y = "") +
  geom_text(aes(label = Description, x = 0, y = Description), hjust = 0, size = 2.5) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(color = "black", size = 7),
    axis.text.y = element_blank(),
    axis.ticks.length.y = unit(0, "in")
  )

ggsave(
  filename = paste0(figure_file_path, 'Figure4/Figure4D.pdf'),
  plot = figure4D,
  height = 1.5, width = 2, units = 'in'
)

# =============================================================================
# Figure 4E - GO Enrichment Analysis for Jurkat
# =============================================================================
# Functional enrichment analysis for up and down regulated O-GlcNAc proteins in Jurkat

# Upregulated O-GlcNAc proteins in Jurkat
OGlcNAc_up_Jurkat <- OGlcNAc_protein_DE_Jurkat |>
  filter(logFC > 0.5, adj.P.Val < 0.05) |>
  pull(Protein.ID)

# Downregulated O-GlcNAc proteins in Jurkat
OGlcNAc_down_Jurkat <- OGlcNAc_protein_DE_Jurkat |>
  filter(logFC < -0.5, adj.P.Val < 0.05) |>
  pull(Protein.ID)

cat("\nJurkat - Up:", length(OGlcNAc_up_Jurkat), ", Down:", length(OGlcNAc_down_Jurkat), "\n")

# GO enrichment for upregulated proteins
OGlcNAc_up_Jurkat_GO <- enrichGO(
  gene = OGlcNAc_up_Jurkat,
  OrgDb = org.Hs.eg.db,
  universe = OGlcNAc_protein_universe,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  OGlcNAc_up_Jurkat_GO@result,
  file = paste0(source_file_path, 'enrichment/OGlcNAc_up_Jurkat_GO.csv')
)

# GO enrichment for downregulated proteins
OGlcNAc_down_Jurkat_GO <- enrichGO(
  gene = OGlcNAc_down_Jurkat,
  OrgDb = org.Hs.eg.db,
  universe = OGlcNAc_protein_universe,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  OGlcNAc_down_Jurkat_GO@result,
  file = paste0(source_file_path, 'enrichment/OGlcNAc_down_Jurkat_GO.csv')
)

cat("Jurkat GO enrichment results saved\n")

# Import saved results for barplot
OGlcNAc_up_Jurkat_GO_result <- read_csv(
  paste0(source_file_path, 'enrichment/OGlcNAc_up_Jurkat_GO.csv')
)
OGlcNAc_down_Jurkat_GO_result <- read_csv(
  paste0(source_file_path, 'enrichment/OGlcNAc_down_Jurkat_GO.csv')
)

# Preview top significant terms
cat("\nJurkat Up - Top 10 significant terms:\n")
OGlcNAc_up_Jurkat_GO_result |>
  filter(pvalue < 0.05) |>
  arrange(pvalue) |>
  dplyr::select(ONTOLOGY, Description, pvalue, Count) |>
  head(10) |>
  print()

cat("\nJurkat Down - Top 10 significant terms:\n")
OGlcNAc_down_Jurkat_GO_result |>
  filter(pvalue < 0.05) |>
  arrange(pvalue) |>
  dplyr::select(ONTOLOGY, Description, pvalue, Count) |>
  head(10) |>
  print()

# Barplot for Jurkat (update Description filter after reviewing results)
# TODO: Update the Description list after manual review of enrichment results
figure4E <- OGlcNAc_down_Jurkat_GO_result |>
  filter(
    Description %in% c(
      'ERAD pathway',
      'proteasomal protein catabolic process',
      'response to endoplasmic reticulum stress',
      'regulation of translational initiation',
      'protein glycosylation'
    )
  ) |>
  mutate(Description = str_replace(Description, "endoplasmic reticulum", "ER")) |>
  ggplot() +
  geom_bar(
    aes(x = -log10(pvalue), y = fct_reorder(Description, -log10(pvalue))),
    stat = "identity", fill = unname(colors_cell["Jurkat"]), width = 0.7
  ) +
  labs(x = expression(-log[10]*"("*paste(italic(P), " Value")*")"), y = "") +
  geom_text(aes(label = Description, x = 0, y = Description), hjust = 0, size = 2.5) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(color = "black", size = 7),
    axis.text.y = element_blank(),
    axis.ticks.length.y = unit(0, "in")
  )

ggsave(
  filename = paste0(figure_file_path, 'Figure4/Figure4E.pdf'),
  plot = figure4E,
  height = 1.5, width = 2, units = 'in'
)

cat("\nFigure 4C-4E: GO enrichment analyses complete.\n")
cat("Results saved to:", source_file_path, "enrichment/\n")
cat("Files: OGlcNAc_up_HEK293T_GO.csv, OGlcNAc_down_HEK293T_GO.csv\n")
cat("       OGlcNAc_up_HepG2_GO.csv, OGlcNAc_down_HepG2_GO.csv\n")
cat("       OGlcNAc_up_Jurkat_GO.csv, OGlcNAc_down_Jurkat_GO.csv\n")
cat("\nReview the top significant terms printed above, then uncomment and update\n")
cat("the barplot code with your selected Description terms.\n")
