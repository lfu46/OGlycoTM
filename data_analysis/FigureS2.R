# Figure S2: Euler plot for localized O-GlcNAc sites across three cell types
# Shows overlap of O-GlcNAc modification sites between HEK293T, HepG2, and Jurkat

library(tidyverse)
library(eulerr)

# Source data paths and colors
source('data_source.R')

# Load O-GlcNAc site data
OGlcNAc_site_HEK293T <- read_csv(
  paste0(source_file_path, 'site/OGlcNAc_site_HEK293T.csv')
)
OGlcNAc_site_HepG2 <- read_csv(
  paste0(source_file_path, 'site/OGlcNAc_site_HepG2.csv')
)
OGlcNAc_site_Jurkat <- read_csv(
  paste0(source_file_path, 'site/OGlcNAc_site_Jurkat.csv')
)

# Extract unique sites from each cell type using site_index
sites_HEK293T <- OGlcNAc_site_HEK293T |>
  pull(site_index) |>
  unique()

sites_HepG2 <- OGlcNAc_site_HepG2 |>
  pull(site_index) |>
  unique()

sites_Jurkat <- OGlcNAc_site_Jurkat |>
  pull(site_index) |>
  unique()

# Print site counts
cat("O-GlcNAc sites per cell type:\n")
cat("HEK293T:", length(sites_HEK293T), "\n")
cat("HepG2:", length(sites_HepG2), "\n")
cat("Jurkat:", length(sites_Jurkat), "\n")

# Create list for euler diagram with cell type names and counts
# Order: HepG2, Jurkat, HEK293T (same as Figure 2A/2B for consistent positioning)
OGlcNAc_site_list <- list(
  sites_HepG2,
  sites_Jurkat,
  sites_HEK293T
)
names(OGlcNAc_site_list) <- c(
  paste0("HepG2 (", length(sites_HepG2), ")"),
  paste0("Jurkat (", length(sites_Jurkat), ")"),
  paste0("HEK293T (", length(sites_HEK293T), ")")
)

# Create euler object (proportional)
OGlcNAc_site_euler <- euler(OGlcNAc_site_list)

# Create color vector (matching the order in the list)
# HepG2 = salmon, Jurkat = green, HEK293T = blue
FigureS2_colors <- c(
  "#F39B7F",
  "#00A087",
  "#4DBBD5"
)
names(FigureS2_colors) <- names(OGlcNAc_site_list)

# Euler plot with labels and quantities
FigureS2 <- plot(
  OGlcNAc_site_euler,
  fills = list(fill = FigureS2_colors, alpha = 0.5),
  edges = list(col = "white", lwd = 2),
  labels = list(font = 1, cex = 0.8),
  quantities = list(font = 1, cex = 0.8)
)

# Save plot
pdf(paste0(figure_file_path, "FigureS2/FigureS2_OGlcNAc_sites_euler.pdf"), width = 2, height = 1.5)
print(FigureS2)
dev.off()

# Display in RStudio
print(FigureS2)

# Print overlap statistics
cat("\nOverlap statistics:\n")
common_sites <- intersect(intersect(sites_HEK293T, sites_HepG2), sites_Jurkat)
cat("Common to all three:", length(common_sites), "\n")

cat("HEK293T & HepG2:", length(intersect(sites_HEK293T, sites_HepG2)), "\n")
cat("HEK293T & Jurkat:", length(intersect(sites_HEK293T, sites_Jurkat)), "\n")
cat("HepG2 & Jurkat:", length(intersect(sites_HepG2, sites_Jurkat)), "\n")

unique_HEK293T <- setdiff(setdiff(sites_HEK293T, sites_HepG2), sites_Jurkat)
unique_HepG2 <- setdiff(setdiff(sites_HepG2, sites_HEK293T), sites_Jurkat)
unique_Jurkat <- setdiff(setdiff(sites_Jurkat, sites_HEK293T), sites_HepG2)
cat("Unique to HEK293T:", length(unique_HEK293T), "\n")
cat("Unique to HepG2:", length(unique_HepG2), "\n")
cat("Unique to Jurkat:", length(unique_Jurkat), "\n")

cat("\nFigure S2 saved to:", figure_file_path, "FigureS2/\n")
