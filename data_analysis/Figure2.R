# import packages
library(tidyverse)

# Figure 2. Identification of O-GlcNAc and O-GalNAc sites and proteins in HEK293T, HepG2 and Jurkat cells


# Figure 2A ---------------------------------------------------------------
# Identification of O-GlcNAc proteins in three types of cells
# Euler plot (proportional Venn diagram) with GO term enrichment plots
# Annotation lines connecting venn regions to GO plots

library(tidyverse)
library(eulerr)
library(cowplot)
library(grid)

# source file path
source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'

# figure file path
figure_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/'

# define color palette
colors_cell <- c("HEK293T" = "#4DBBD5", "HepG2" = "#F39B7F", "Jurkat" = "#00A087")

# load bonafide glyco data
OGlyco_HEK293T_bonafide <- read_csv(
  paste0(source_file_path, 'filtered/OGlyco_HEK293T_bonafide.csv')
)
OGlyco_HepG2_bonafide <- read_csv(
  paste0(source_file_path, 'filtered/OGlyco_HepG2_bonafide.csv')
)
OGlyco_Jurkat_bonafide <- read_csv(
  paste0(source_file_path, 'filtered/OGlyco_Jurkat_bonafide.csv')
)

# Extract O-GlcNAc proteins from each cell type
OGlcNAc_proteins_HEK293T <- OGlyco_HEK293T_bonafide |>
  filter(Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859')) |>
  pull(Protein.ID) |>
  unique()

OGlcNAc_proteins_HepG2 <- OGlyco_HepG2_bonafide |>
  filter(Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859')) |>
  pull(Protein.ID) |>
  unique()

OGlcNAc_proteins_Jurkat <- OGlyco_Jurkat_bonafide |>
  filter(Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859')) |>
  pull(Protein.ID) |>
  unique()

# Create list for Euler diagram (no cell type labels, just for calculation)
# Order: HepG2, Jurkat, HEK293T to position them correctly
# HepG2 (left), HEK293T (upper-right), Jurkat (lower-right)
OGlcNAc_protein_list <- list(
  HepG2 = OGlcNAc_proteins_HepG2,
  Jurkat = OGlcNAc_proteins_Jurkat,
  HEK293T = OGlcNAc_proteins_HEK293T
)

# Create euler object (proportional)
OGlcNAc_protein_euler <- euler(OGlcNAc_protein_list)

# Create color vector (matching the order in the list)
# HepG2 = salmon (#F39B7F), Jurkat = green (#00A087), HEK293T = blue (#4DBBD5)
Figure2A_colors <- c(
  "HepG2" = "#F39B7F",
  "Jurkat" = "#00A087",
  "HEK293T" = "#4DBBD5"
)

# Euler diagram (no labels, only quantities/numbers shown)
Figure2A_euler <- plot(
  OGlcNAc_protein_euler,
  fills = list(fill = Figure2A_colors, alpha = 0.5),
  edges = list(col = "white", lwd = 2),
  labels = FALSE,
  quantities = list(font = 1, cex = 1.0)
)

# Calculate protein counts for each category
# Common (intersection of all 3)
common_OGlcNAc_proteins <- intersect(
  intersect(OGlcNAc_proteins_HEK293T, OGlcNAc_proteins_HepG2),
  OGlcNAc_proteins_Jurkat
)

# Unique to each cell type
unique_OGlcNAc_HEK293T <- setdiff(
  setdiff(OGlcNAc_proteins_HEK293T, OGlcNAc_proteins_HepG2),
  OGlcNAc_proteins_Jurkat
)
unique_OGlcNAc_HepG2 <- setdiff(
  setdiff(OGlcNAc_proteins_HepG2, OGlcNAc_proteins_HEK293T),
  OGlcNAc_proteins_Jurkat
)
unique_OGlcNAc_Jurkat <- setdiff(
  setdiff(OGlcNAc_proteins_Jurkat, OGlcNAc_proteins_HEK293T),
  OGlcNAc_proteins_HepG2
)

# Total O-GlcNAc proteins (for enrichment universe)
total_OGlcNAc_proteins <- union(
  union(OGlcNAc_proteins_HEK293T, OGlcNAc_proteins_HepG2),
  OGlcNAc_proteins_Jurkat
)

# --- Load GO enrichment results ---
# Common O-GlcNAc GO terms (select 6 key terms)
common_GO_selected <- read_csv(
  paste0(source_file_path, 'enrichment/Figure2C_common_OGlcNAc_GO.csv')
) |>
  filter(
    Description %in% c(
      'transcription coregulator activity',
      'mRNA binding',
      'nucleocytoplasmic transport',
      'nuclear transport',
      'cytoplasmic stress granule',
      'histone acetyltransferase complex'
    )
  ) |>
  dplyr::select(Description, pvalue, Count) |>
  mutate(
    log_pvalue = -log10(pvalue),
    # Abbreviate long terms
    Description = case_when(
      Description == "transcription coregulator activity" ~ "Transcription coregulator",
      Description == "nucleocytoplasmic transport" ~ "Nucleocytoplasmic transport",
      Description == "histone acetyltransferase complex" ~ "Histone acetyltransferase",
      Description == "cytoplasmic stress granule" ~ "Cytoplasmic stress granule",
      TRUE ~ Description
    )
  )

# HEK293T unique GO terms
HEK293T_GO_selected <- read_csv(
  paste0(source_file_path, 'enrichment/Figure2C_unique_HEK293T_OGlcNAc_GO.csv')
) |>
  filter(
    Description %in% c(
      'chaperone-mediated protein folding',
      'cell cycle process',
      'epithelial cell proliferation'
    )
  ) |>
  dplyr::select(Description, pvalue, Count) |>
  mutate(
    log_pvalue = -log10(pvalue),
    # Abbreviate long terms
    Description = case_when(
      Description == "chaperone-mediated protein folding" ~ "Chaperone-mediated folding",
      Description == "epithelial cell proliferation" ~ "Epithelial cell proliferation",
      TRUE ~ Description
    )
  )

# HepG2 unique GO terms
HepG2_GO_selected <- read_csv(
  paste0(source_file_path, 'enrichment/Figure2C_unique_HepG2_OGlcNAc_GO.csv')
) |>
  filter(
    Description %in% c(
      'extracellular matrix',
      'cell periphery',
      'endoplasmic reticulum lumen'
    )
  ) |>
  dplyr::select(Description, pvalue, Count) |>
  mutate(
    log_pvalue = -log10(pvalue),
    # Abbreviate long terms
    Description = case_when(
      Description == "endoplasmic reticulum lumen" ~ "ER lumen",
      Description == "extracellular matrix" ~ "Extracellular matrix",
      TRUE ~ Description
    )
  )

# Jurkat unique GO terms
Jurkat_GO_selected <- read_csv(
  paste0(source_file_path, 'enrichment/Figure2C_unique_Jurkat_OGlcNAc_GO.csv')
) |>
  filter(
    Description %in% c(
      'leukocyte activation involved in immune response',
      'leukocyte cell-cell adhesion',
      'T cell activation'
    )
  ) |>
  dplyr::select(Description, pvalue, Count) |>
  mutate(
    log_pvalue = -log10(pvalue),
    # Abbreviate long terms
    Description = case_when(
      Description == "leukocyte activation involved in immune response" ~ "Leukocyte activation",
      Description == "leukocyte cell-cell adhesion" ~ "Leukocyte adhesion",
      TRUE ~ Description
    )
  )

# --- Create individual GO term bar plots with gradient effect ---
# Function to create gradient bar data
create_gradient_data <- function(df, n_segments = 50) {
  df <- df |> mutate(y_num = as.numeric(fct_reorder(Description, log_pvalue)))

  gradient_df <- do.call(rbind, lapply(1:nrow(df), function(i) {
    row <- df[i, ]
    segments <- data.frame(
      Description = row$Description,
      y_num = row$y_num,
      xmin = seq(0, row$log_pvalue, length.out = n_segments + 1)[-(n_segments + 1)],
      xmax = seq(0, row$log_pvalue, length.out = n_segments + 1)[-1],
      segment = 1:n_segments
    )
    segments$alpha_val <- seq(0.9, 0.3, length.out = n_segments)
    segments
  }))
  gradient_df
}

# Lighter pastel colors for gradient
color_common <- "#A0A0A0"
color_HEK293T <- "#7DCDE5"
color_HepG2 <- "#F7BBA8"
color_Jurkat <- "#4DC4B0"

# Common GO plot (top left) with gradient
common_gradient_df <- create_gradient_data(common_GO_selected)
common_GO_plot <- ggplot() +
  geom_rect(
    data = common_gradient_df,
    aes(xmin = xmin, xmax = xmax,
        ymin = y_num - 0.35, ymax = y_num + 0.35,
        alpha = alpha_val),
    fill = color_common
  ) +
  geom_text(
    data = common_GO_selected |> mutate(y_num = as.numeric(fct_reorder(Description, log_pvalue))),
    aes(label = Description, x = 0.05, y = y_num),
    hjust = 0, size = 2.8, color = "black"
  ) +
  scale_alpha_identity() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_continuous(breaks = 1:nrow(common_GO_selected), labels = NULL) +
  labs(x = expression(-log[10]~"("*italic(p)~value*")"), y = NULL) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica", color = "black"),
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(color = "black", size = 7),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank()
  )

# HEK293T GO plot (top right) with gradient
HEK293T_gradient_df <- create_gradient_data(HEK293T_GO_selected)
HEK293T_GO_plot <- ggplot() +
  geom_rect(
    data = HEK293T_gradient_df,
    aes(xmin = xmin, xmax = xmax,
        ymin = y_num - 0.35, ymax = y_num + 0.35,
        alpha = alpha_val),
    fill = color_HEK293T
  ) +
  geom_text(
    data = HEK293T_GO_selected |> mutate(y_num = as.numeric(fct_reorder(Description, log_pvalue))),
    aes(label = Description, x = 0.05, y = y_num),
    hjust = 0, size = 2.8, color = "black"
  ) +
  scale_alpha_identity() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_continuous(breaks = 1:nrow(HEK293T_GO_selected), labels = NULL) +
  labs(x = expression(-log[10]~"("*italic(p)~value*")"), y = NULL) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica", color = "black"),
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(color = "black", size = 7),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank()
  )

# HepG2 GO plot (bottom left) with gradient
HepG2_gradient_df <- create_gradient_data(HepG2_GO_selected)
HepG2_GO_plot <- ggplot() +
  geom_rect(
    data = HepG2_gradient_df,
    aes(xmin = xmin, xmax = xmax,
        ymin = y_num - 0.35, ymax = y_num + 0.35,
        alpha = alpha_val),
    fill = color_HepG2
  ) +
  geom_text(
    data = HepG2_GO_selected |> mutate(y_num = as.numeric(fct_reorder(Description, log_pvalue))),
    aes(label = Description, x = 0.05, y = y_num),
    hjust = 0, size = 2.8, color = "black"
  ) +
  scale_alpha_identity() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_continuous(breaks = 1:nrow(HepG2_GO_selected), labels = NULL) +
  labs(x = expression(-log[10]~"("*italic(p)~value*")"), y = NULL) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica", color = "black"),
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(color = "black", size = 7),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank()
  )

# Jurkat GO plot (bottom right) with gradient
Jurkat_gradient_df <- create_gradient_data(Jurkat_GO_selected)
Jurkat_GO_plot <- ggplot() +
  geom_rect(
    data = Jurkat_gradient_df,
    aes(xmin = xmin, xmax = xmax,
        ymin = y_num - 0.35, ymax = y_num + 0.35,
        alpha = alpha_val),
    fill = color_Jurkat
  ) +
  geom_text(
    data = Jurkat_GO_selected |> mutate(y_num = as.numeric(fct_reorder(Description, log_pvalue))),
    aes(label = Description, x = 0.05, y = y_num),
    hjust = 0, size = 2.8, color = "black"
  ) +
  scale_alpha_identity() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_continuous(breaks = 1:nrow(Jurkat_GO_selected), labels = NULL) +
  labs(x = expression(-log[10]~"("*italic(p)~value*")"), y = NULL) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica", color = "black"),
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(color = "black", size = 7),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.background = element_blank(),
    panel.background = element_blank()
  )

# --- Create the combined figure using cowplot ---
# Convert euler plot to grob for cowplot
euler_grob <- grid::grid.grabExpr(print(Figure2A_euler), width = 3, height = 3)

# Create the combined layout
# Layout: Left column (Jurkat, HEK293T) | Center (Euler, smaller) | Right column (HepG2, Common)
# Labels in the middle area between venn and GO plots

Figure2A <- ggdraw() +
  # Central euler diagram
  draw_grob(euler_grob, x = 0.30, y = 0.20, width = 0.40, height = 0.60) +
  # Left column - Common (larger, 6 terms) and HepG2 (smaller, 3 terms)
  # Common GO plot (top left) - larger height for 6 terms
  draw_plot(common_GO_plot, x = 0, y = 0.48, width = 0.30, height = 0.52) +
  # HepG2 GO plot (bottom left) - smaller height for 3 terms
  draw_plot(HepG2_GO_plot, x = 0, y = 0.12, width = 0.30, height = 0.36) +
  # Right column - HEK293T and Jurkat (smaller, 3 terms each)
  # HEK293T GO plot (top right) - smaller height, centered vertically
  draw_plot(HEK293T_GO_plot, x = 0.70, y = 0.58, width = 0.30, height = 0.36) +
  # Jurkat GO plot (bottom right) - smaller height, centered vertically
  draw_plot(Jurkat_GO_plot, x = 0.70, y = 0.06, width = 0.30, height = 0.36) +
  # Labels in middle area (between venn and GO plots) with protein counts
  # Use darker/more saturated colors for labels
  # Common label (top left connector area)
  draw_label(paste0("Common (", length(common_OGlcNAc_proteins), ")"),
             x = 0.35, y = 0.94, hjust = 0.5, vjust = 0.5,
             fontface = "plain", color = "#505050", size = 7) +
  # HepG2 label (bottom left connector area)
  draw_label(paste0("HepG2 (", length(unique_OGlcNAc_HepG2), ")"),
             x = 0.35, y = 0.06, hjust = 0.5, vjust = 0.5,
             fontface = "plain", color = "#D97A5A", size = 7) +
  # HEK293T label (top right connector area)
  draw_label(paste0("HEK293T (", length(unique_OGlcNAc_HEK293T), ")"),
             x = 0.65, y = 0.94, hjust = 0.5, vjust = 0.5,
             fontface = "plain", color = "#2A9BB5", size = 7) +
  # Jurkat label (bottom right connector area)
  draw_label(paste0("Jurkat (", length(unique_OGlcNAc_Jurkat), ")"),
             x = 0.65, y = 0.06, hjust = 0.5, vjust = 0.5,
             fontface = "plain", color = "#008570", size = 7) +
  # Annotation lines with circles (using darker colors to match labels)
  # Line from Common (402) to Common label (top left)
  geom_path(
    data = data.frame(x = c(0.50, 0.42, 0.35), y = c(0.50, 0.82, 0.90)),
    aes(x = x, y = y), color = "#505050", linewidth = 0.5
  ) +
  # Circle around Common count
  geom_point(aes(x = 0.50, y = 0.50), shape = 1, size = 5, color = "#505050", stroke = 0.8) +
  # Line from HepG2 unique (145) to HepG2 label (bottom left)
  geom_path(
    data = data.frame(x = c(0.42, 0.38, 0.35), y = c(0.38, 0.16, 0.10)),
    aes(x = x, y = y), color = "#D97A5A", linewidth = 0.5
  ) +
  # Circle around HepG2 count
  geom_point(aes(x = 0.42, y = 0.38), shape = 1, size = 4, color = "#D97A5A", stroke = 0.8) +
  # Line from HEK293T unique (187) to HEK293T label (top right)
  geom_path(
    data = data.frame(x = c(0.58, 0.62, 0.65), y = c(0.68, 0.84, 0.90)),
    aes(x = x, y = y), color = "#2A9BB5", linewidth = 0.5
  ) +
  # Circle around HEK293T count
  geom_point(aes(x = 0.58, y = 0.68), shape = 1, size = 4, color = "#2A9BB5", stroke = 0.8) +
  # Line from Jurkat unique (145) to Jurkat label (bottom right)
  geom_path(
    data = data.frame(x = c(0.58, 0.62, 0.65), y = c(0.32, 0.16, 0.10)),
    aes(x = x, y = y), color = "#008570", linewidth = 0.5
  ) +
  # Circle around Jurkat count
  geom_point(aes(x = 0.58, y = 0.32), shape = 1, size = 4, color = "#008570", stroke = 0.8)

# Save Figure 2A
ggsave(
  filename = paste0(figure_file_path, "Figure2/Figure2A.pdf"),
  plot = Figure2A,
  width = 6, height = 3, units = "in"
)


# Figure 2B ---------------------------------------------------------------
# Euler plot for O-GalNAc proteins across three cell types

library(tidyverse)
library(eulerr)

# source file path
source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'

# figure file path
figure_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/'

# load bonafide glyco data
OGlyco_HEK293T_bonafide <- read_csv(
  paste0(source_file_path, 'filtered/OGlyco_HEK293T_bonafide.csv')
)
OGlyco_HepG2_bonafide <- read_csv(
  paste0(source_file_path, 'filtered/OGlyco_HepG2_bonafide.csv')
)
OGlyco_Jurkat_bonafide <- read_csv(
  paste0(source_file_path, 'filtered/OGlyco_Jurkat_bonafide.csv')
)

# Extract O-GalNAc proteins from each cell type
# O-GalNAc criteria: GAO_Methoxylamine modification with mass shifts 326.1339 and 555.2968
OGalNAc_proteins_HEK293T <- OGlyco_HEK293T_bonafide |>
  filter(Total.Glycan.Composition %in% c(
    'HexNAt(1)GAO_Methoxylamine(1) % 326.1339',
    'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968'
  )) |>
  pull(Protein.ID) |>
  unique()

OGalNAc_proteins_HepG2 <- OGlyco_HepG2_bonafide |>
  filter(Total.Glycan.Composition %in% c(
    'HexNAt(1)GAO_Methoxylamine(1) % 326.1339',
    'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968'
  )) |>
  pull(Protein.ID) |>
  unique()

OGalNAc_proteins_Jurkat <- OGlyco_Jurkat_bonafide |>
  filter(Total.Glycan.Composition %in% c(
    'HexNAt(1)GAO_Methoxylamine(1) % 326.1339',
    'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968'
  )) |>
  pull(Protein.ID) |>
  unique()

# Create list for euler diagram with cell type names and counts
# Order: HepG2, Jurkat, HEK293T (same as Figure 2A for consistent positioning)
OGalNAc_protein_list <- list(
  OGalNAc_proteins_HepG2,
  OGalNAc_proteins_Jurkat,
  OGalNAc_proteins_HEK293T
)
names(OGalNAc_protein_list) <- c(

  paste0("HepG2 (", length(OGalNAc_proteins_HepG2), ")"),
  paste0("Jurkat (", length(OGalNAc_proteins_Jurkat), ")"),
  paste0("HEK293T (", length(OGalNAc_proteins_HEK293T), ")")
)

# Create euler object (proportional)
OGalNAc_protein_euler <- euler(OGalNAc_protein_list)

# Create color vector (matching the order in the list)
# HepG2 = salmon, Jurkat = green, HEK293T = blue
Figure2B_colors <- c(
  "#F39B7F",
  "#00A087",
  "#4DBBD5"
)
names(Figure2B_colors) <- names(OGalNAc_protein_list)

# Euler plot with labels and quantities
Figure2B <- plot(
  OGalNAc_protein_euler,
  fills = list(fill = Figure2B_colors, alpha = 0.5),
  edges = list(col = "white", lwd = 2),
  labels = list(font = 1, cex = 0.8),
  quantities = list(font = 1, cex = 0.8)
)

# Save plot
pdf(paste0(figure_file_path, "Figure2/Figure2B.pdf"), width = 2, height = 1.5)
print(Figure2B)
dev.off()


# Figure 2C - Enrichment ---------------------------------------------------
# GO enrichment for exclusive O-GalNAc proteins
# Run this section first to generate enrichment results

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

# source file path
source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'

# load bonafide glyco data
OGlyco_HEK293T_bonafide <- read_csv(
  paste0(source_file_path, 'filtered/OGlyco_HEK293T_bonafide.csv')
)
OGlyco_HepG2_bonafide <- read_csv(
  paste0(source_file_path, 'filtered/OGlyco_HepG2_bonafide.csv')
)
OGlyco_Jurkat_bonafide <- read_csv(
  paste0(source_file_path, 'filtered/OGlyco_Jurkat_bonafide.csv')
)

# Define glycan composition filters
OGlcNAc_compositions <- c(
  'HexNAt(1) % 299.1230',
  'HexNAt(1)TMT6plex(1) % 528.2859'
)
OGalNAc_compositions <- c(
  'HexNAt(1)GAO_Methoxylamine(1) % 326.1339',
  'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968'
)

# Combine all cell types
all_bonafide <- bind_rows(
  OGlyco_HEK293T_bonafide,
  OGlyco_HepG2_bonafide,
  OGlyco_Jurkat_bonafide
)

# Identify O-GlcNAc and O-GalNAc proteins (across all cell types)
OGlcNAc_proteins <- all_bonafide |>
  filter(Total.Glycan.Composition %in% OGlcNAc_compositions) |>
  pull(Protein.ID) |>
  unique()

OGalNAc_proteins <- all_bonafide |>
  filter(Total.Glycan.Composition %in% OGalNAc_compositions) |>
  pull(Protein.ID) |>
  unique()

# Identify exclusive proteins
exclusive_OGalNAc_proteins <- setdiff(OGalNAc_proteins, OGlcNAc_proteins)
exclusive_OGlcNAc_proteins <- setdiff(OGlcNAc_proteins, OGalNAc_proteins)

# Create universe: exclusive O-GlcNAc + exclusive O-GalNAc proteins only
# (excluding dual-labeled proteins)
universe_exclusive <- union(exclusive_OGlcNAc_proteins, exclusive_OGalNAc_proteins)

cat("Exclusive O-GlcNAc proteins:", length(exclusive_OGlcNAc_proteins), "\n")
cat("Exclusive O-GalNAc proteins:", length(exclusive_OGalNAc_proteins), "\n")
cat("Universe (exclusive glycoproteins only):", length(universe_exclusive), "\n")

# GO enrichment for exclusive O-GalNAc proteins
exclusive_OGalNAc_GO <- enrichGO(
  gene = exclusive_OGalNAc_proteins,
  OrgDb = org.Hs.eg.db,
  universe = universe_exclusive,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

# Save enrichment results to CSV
write_csv(
  exclusive_OGalNAc_GO@result,
  file = paste0(source_file_path, 'enrichment/exclusive_OGalNAc_GO.csv')
)

cat("\nGO enrichment results saved to:", source_file_path, "enrichment/exclusive_OGalNAc_GO.csv\n")

# Preview top significant terms
cat("\nTop 20 significant GO terms (p < 0.05):\n")
exclusive_OGalNAc_GO@result |>
  filter(pvalue < 0.05) |>
  arrange(pvalue) |>
  dplyr::select(ONTOLOGY, Description, pvalue, Count) |>
  head(20) |>
  print()


# Figure 2C - Dotplot ------------------------------------------------------
# Create dotplot with manually selected GO terms
# Uses gradient-filled circles for visual effect

library(tidyverse)
library(ggforce)

# source file path
source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'

# figure file path
figure_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/'

# Import saved enrichment results
exclusive_OGalNAc_GO_result <- read_csv(
  paste0(source_file_path, 'enrichment/exclusive_OGalNAc_GO.csv')
)

# Manually select 5-6 GO terms (fill in after reviewing results)
Figure2C_GO_selected <- exclusive_OGalNAc_GO_result |>
  filter(
    Description %in% c(
      # TODO: Add 5-6 GO terms here after reviewing enrichment results
      'cell surface',
      'membrane',
      'extracellular region',
      'transporter activity',
      'glycosaminoglycan binding',
      'Golgi apparatus'
    )
  ) |>
  dplyr::select(ONTOLOGY, Description, pvalue, Count) |>
  mutate(
    log_pvalue = -log10(pvalue),
    Description = str_wrap(Description, width = 40)
  )

# Define color for O-GalNAc
color_OGalNAc <- "#4DBBD5"

# Prepare data with numeric y-axis and scaled radius
Figure2C_plot_data <- Figure2C_GO_selected |>
  mutate(
    y_num = as.numeric(fct_reorder(Description, log_pvalue)),
    # Scale radius based on Count (larger range for bigger dots)
    radius_y = scales::rescale(Count, to = c(0.25, 0.50))
  )

# Calculate aspect ratio correction factor
# Based on data range and plot dimensions (width=3, height=2)
x_range <- max(Figure2C_plot_data$log_pvalue) - 0
y_range <- max(Figure2C_plot_data$y_num) - min(Figure2C_plot_data$y_num) + 2  # add padding
plot_ratio <- 0.8 * (3 / 2) * (y_range / x_range)  # adjusted for circular appearance

# Add x-radius adjusted for aspect ratio
Figure2C_plot_data <- Figure2C_plot_data |>
  mutate(radius_x = radius_y * plot_ratio)

# Create gradient ellipse data - concentric rings with color gradient
# Gradient goes from white (center) to main color (edge)
n_rings <- 20
gradient_ellipse_data <- Figure2C_plot_data |>
  rowwise() |>
  mutate(
    rings = list(tibble(
      ring = 1:n_rings,
      a = seq(radius_x, radius_x * 0.05, length.out = n_rings),  # x-radius
      b = seq(radius_y, radius_y * 0.05, length.out = n_rings),  # y-radius
      # Color interpolation from main color (outer) to white (inner)
      ring_color = colorRampPalette(c(color_OGalNAc, "white"))(n_rings)
    ))
  ) |>
  unnest(rings)

# Create dotplot with gradient ellipses (appear as circles due to aspect ratio correction)
Figure2C <- ggplot() +
  # Dashed lines connecting to y-axis
  geom_segment(
    data = Figure2C_plot_data,
    aes(x = 0, xend = log_pvalue, y = y_num, yend = y_num),
    linetype = "dashed",
    color = "grey70",
    linewidth = 0.3
  ) +
  # Gradient ellipses (draw outer rings first, inner rings on top)
  geom_ellipse(
    data = gradient_ellipse_data,
    aes(x0 = log_pvalue, y0 = y_num, a = a, b = b, angle = 0, fill = ring_color),
    color = NA
  ) +
  # Invisible points just for the size legend
  geom_point(
    data = Figure2C_plot_data,
    aes(x = log_pvalue, y = y_num, size = Count),
    alpha = 0
  ) +
  scale_fill_identity() +
  scale_size_continuous(range = c(2, 6), breaks = c(10, 20, 30, 40)) +
  # Custom y-axis with original descriptions
  scale_y_continuous(
    breaks = Figure2C_plot_data$y_num,
    labels = Figure2C_plot_data$Description,
    expand = expansion(mult = c(0.15, 0.15))
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = expression(-log[10]~"("*italic(p)~value*")"),
    y = NULL,
    size = "Count"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica", color = "black"),
    axis.text = element_text(color = "black", size = 9),
    axis.text.y = element_text(lineheight = 0.8),
    axis.title = element_text(size = 9),
    legend.position = "right",
    legend.key.size = unit(0.4, 'cm'),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  ) +
  guides(size = guide_legend(override.aes = list(alpha = 1, fill = color_OGalNAc)))

# Save plot
ggsave(
  filename = paste0(figure_file_path, "Figure2/Figure2C.pdf"),
  plot = Figure2C,
  width = 3, height = 2, units = "in"
)

cat("\nFigure 2C saved to:", figure_file_path, "Figure2/Figure2C.pdf\n")


# Figure 2D ---------------------------------------------------------------
# Statistical summary of 144/138 ratio for ALL glycoproteins by subcellular location
# Uses PSM level data from all three cell types combined
# Proteins categorized by HPA annotation (unique annotation only)

library(tidyverse)
library(ggpubr)

# source file path
source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'

# figure file path
figure_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/'

cat("=== Figure 2D: 144/138 Ratio by Subcellular Location (All Glycoproteins) ===\n\n")

# Load Human Protein Atlas subcellular location data
subcellular_location <- read_tsv(
  paste0(source_file_path, 'reference/subcellular_location.tsv'),
  show_col_types = FALSE
)

# Filter for single location and reliable annotation (unique annotation only)
subcellular_location_filtered <- subcellular_location |>
  filter(!str_detect(`Main location`, ";")) |>
  filter(Reliability != "Uncertain") |>
  dplyr::select(Gene_name = `Gene name`, Location = `Main location`, Reliability) |>
  group_by(Gene_name) |>
  filter(n() == 1) |>
  ungroup() |>
  distinct(Gene_name, Location, .keep_all = TRUE)

cat("Filtered subcellular location data:", nrow(subcellular_location_filtered), "genes\n")

# Map HPA locations to reference categories
location_mapping <- c(
  "Nucleoplasm" = "Nucleus",
  "Nucleoli" = "Nucleus",
  "Nuclear membrane" = "Nucleus",
  "Nuclear bodies" = "Nucleus",
  "Nuclear speckles" = "Nucleus",
  "Cytosol" = "Cytoplasm",
  "Mitochondria" = "Mitochondrion",
  "Endoplasmic reticulum" = "Endoplasmic reticulum",
  "Golgi apparatus" = "Golgi apparatus",
  "Lysosomes" = "Lysosome",
  "Plasma membrane" = "Plasma membrane",
  "Cell Junctions" = "Extracellular"
)

# Add mapped location column
subcellular_location_mapped <- subcellular_location_filtered |>
  filter(Location %in% names(location_mapping)) |>
  mutate(Location_mapped = location_mapping[Location])

cat("\nMapped locations in HPA:\n")
subcellular_location_mapped |>
  group_by(Location_mapped) |>
  summarise(n = n(), .groups = "drop") |>
  arrange(desc(n)) |>
  print()

# Load bonafide glyco data
OGlyco_HEK293T_bonafide <- read_csv(
  paste0(source_file_path, 'filtered/OGlyco_HEK293T_bonafide.csv'),
  show_col_types = FALSE
)
OGlyco_HepG2_bonafide <- read_csv(
  paste0(source_file_path, 'filtered/OGlyco_HepG2_bonafide.csv'),
  show_col_types = FALSE
)
OGlyco_Jurkat_bonafide <- read_csv(
  paste0(source_file_path, 'filtered/OGlyco_Jurkat_bonafide.csv'),
  show_col_types = FALSE
)

# Combine all cell types
all_bonafide <- bind_rows(
  OGlyco_HEK293T_bonafide |> mutate(Cell_Type = "HEK293T"),
  OGlyco_HepG2_bonafide |> mutate(Cell_Type = "HepG2"),
  OGlyco_Jurkat_bonafide |> mutate(Cell_Type = "Jurkat")
)

cat("\nTotal PSMs:", nrow(all_bonafide), "\n")
cat("Unique proteins:", n_distinct(all_bonafide$Protein.ID), "\n")

# Check raw Ratio.138.144 distribution
cat("\n--- Diagnostic: Ratio.138.144 distribution ---\n")
ratio_raw <- as.numeric(all_bonafide$Ratio.138.144)
cat("Total values:", length(ratio_raw), "\n")
cat("NA values:", sum(is.na(ratio_raw)), "\n")
cat("Zero values:", sum(ratio_raw == 0, na.rm = TRUE), "\n")
cat("Negative values:", sum(ratio_raw < 0, na.rm = TRUE), "\n")
cat("Positive values:", sum(ratio_raw > 0, na.rm = TRUE), "\n")
cat("Summary of positive values:\n")
print(summary(ratio_raw[ratio_raw > 0]))

# Extract ALL PSMs with unique HPA location annotation
# Filter out zero, negative, and NA Ratio.138.144 values
all_PSMs_with_location <- all_bonafide |>
  inner_join(
    subcellular_location_mapped |> dplyr::select(Gene_name, Location_mapped),
    by = c("Gene" = "Gene_name")
  ) |>
  mutate(Ratio_138_144_numeric = as.numeric(Ratio.138.144)) |>
  filter(!is.na(Ratio_138_144_numeric) & Ratio_138_144_numeric > 0) |>
  mutate(Ratio_144_138 = 1 / Ratio_138_144_numeric) |>
  filter(is.finite(Ratio_144_138) & Ratio_144_138 > 0)

cat("\n--- After filtering ---\n")
cat("Ratio_144_138 summary:\n")
print(summary(all_PSMs_with_location$Ratio_144_138))
cat("Zero values in Ratio_144_138:", sum(all_PSMs_with_location$Ratio_144_138 == 0), "\n")

cat("\nPSMs with unique HPA location annotation:", nrow(all_PSMs_with_location), "\n")
cat("Proteins with unique HPA location annotation:", n_distinct(all_PSMs_with_location$Protein.ID), "\n")

# Statistical Summary by Subcellular Location
cat("\n", strrep("=", 70), "\n")
cat("STATISTICAL SUMMARY: All Glycoproteins by Subcellular Location\n")
cat(strrep("=", 70), "\n\n")

location_summary <- all_PSMs_with_location |>
  group_by(Location_mapped) |>
  summarise(
    n_PSMs = n(),
    n_proteins = n_distinct(Protein.ID),
    median_ratio = round(median(Ratio_144_138), 4),
    mean_ratio = round(mean(Ratio_144_138), 4),
    sd_ratio = round(sd(Ratio_144_138), 4),
    Q1 = round(quantile(Ratio_144_138, 0.25), 4),
    Q3 = round(quantile(Ratio_144_138, 0.75), 4),
    .groups = "drop"
  ) |>
  arrange(desc(n_PSMs))

print(location_summary, n = Inf)

# Statistical Summary by Glycan Type and Location
OGlcNAc_compositions <- c(
  'HexNAt(1) % 299.1230',
  'HexNAt(1)TMT6plex(1) % 528.2859'
)
OGalNAc_compositions <- c(
  'HexNAt(1)GAO_Methoxylamine(1) % 326.1339',
  'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968'
)

# Add glycan type classification
all_PSMs_with_location <- all_PSMs_with_location |>
  mutate(
    Glycan_Type = case_when(
      Total.Glycan.Composition %in% OGlcNAc_compositions ~ "O-GlcNAc",
      Total.Glycan.Composition %in% OGalNAc_compositions ~ "O-GalNAc",
      TRUE ~ "Other"
    )
  )

cat("\n", strrep("=", 70), "\n")
cat("STATISTICAL SUMMARY: By Glycan Type\n")
cat(strrep("=", 70), "\n\n")

glycan_summary <- all_PSMs_with_location |>
  group_by(Glycan_Type) |>
  summarise(
    n_PSMs = n(),
    n_proteins = n_distinct(Protein.ID),
    median_ratio = round(median(Ratio_144_138), 4),
    mean_ratio = round(mean(Ratio_144_138), 4),
    sd_ratio = round(sd(Ratio_144_138), 4),
    .groups = "drop"
  )

print(glycan_summary)

cat("\n", strrep("=", 70), "\n")
cat("STATISTICAL SUMMARY: By Glycan Type AND Location\n")
cat(strrep("=", 70), "\n\n")

glycan_location_summary <- all_PSMs_with_location |>
  filter(Glycan_Type %in% c("O-GlcNAc", "O-GalNAc")) |>
  group_by(Glycan_Type, Location_mapped) |>
  summarise(
    n_PSMs = n(),
    n_proteins = n_distinct(Protein.ID),
    median_ratio = round(median(Ratio_144_138), 4),
    mean_ratio = round(mean(Ratio_144_138), 4),
    sd_ratio = round(sd(Ratio_144_138), 4),
    .groups = "drop"
  ) |>
  arrange(Glycan_Type, desc(n_PSMs))

print(glycan_location_summary, n = Inf)

# Wide format for easy comparison
cat("\n", strrep("=", 70), "\n")
cat("COMPARISON: O-GlcNAc vs O-GalNAc by Location (median ratio)\n")
cat(strrep("=", 70), "\n\n")

comparison_wide <- glycan_location_summary |>
  dplyr::select(Glycan_Type, Location_mapped, n_PSMs, median_ratio) |>
  pivot_wider(
    names_from = Glycan_Type,
    values_from = c(n_PSMs, median_ratio),
    values_fill = list(n_PSMs = 0, median_ratio = NA)
  ) |>
  arrange(desc(`n_PSMs_O-GlcNAc`))

print(comparison_wide, n = Inf)

# Prepare data for plotting
Figure2D_data <- all_PSMs_with_location |>
  mutate(
    Location_mapped = factor(Location_mapped, levels = c(
      "Nucleus", "Cytoplasm", "Mitochondrion", "Endoplasmic reticulum",
      "Golgi apparatus", "Lysosome", "Plasma membrane", "Extracellular"
    ))
  )

cat("\n", strrep("=", 70), "\n")
cat("SAMPLE SIZES FOR PLOT\n")
cat(strrep("=", 70), "\n\n")

sample_sizes <- Figure2D_data |>
  group_by(Location_mapped) |>
  summarise(n_PSMs = n(), n_proteins = n_distinct(Protein.ID), .groups = "drop") |>
  arrange(match(Location_mapped, c(
    "Nucleus", "Cytoplasm", "Mitochondrion", "Endoplasmic reticulum",
    "Golgi apparatus", "Lysosome", "Plasma membrane", "Extracellular"
  )))

print(sample_sizes)

# Combine Plasma membrane and Extracellular, use abbreviations for plot
Figure2D_plot_data <- Figure2D_data |>
  mutate(
    Location_plot = case_when(
      Location_mapped == "Nucleus" ~ "Nucleus",
      Location_mapped == "Cytoplasm" ~ "Cytoplasm",
      Location_mapped == "Mitochondrion" ~ "Mito",
      Location_mapped == "Endoplasmic reticulum" ~ "ER",
      Location_mapped == "Golgi apparatus" ~ "Golgi",
      Location_mapped == "Plasma membrane" ~ "PM",
      Location_mapped == "Extracellular" ~ "EC",
      TRUE ~ NA_character_
    ),
    Location_plot = factor(Location_plot, levels = c(
      "Nucleus", "Cytoplasm", "Mito", "ER", "Golgi", "PM", "EC"
    ))
  ) |>
  filter(!is.na(Location_plot))

# Calculate sample sizes for annotation
sample_sizes <- Figure2D_plot_data |>
  group_by(Location_plot) |>
  summarise(n = n(), .groups = "drop")

cat("\nSample sizes for plot:\n")
print(sample_sizes)

# Define colors for locations (same as Figure5.R)
colors_location <- c(
  "Nucleus" = "#3C5488",
  "Cytoplasm" = "#00A087",
  "Mito" = "#E64B35",
  "ER" = "#7E6148",
  "Golgi" = "#8491B4",
  "PM" = "#91D1C2",
  "EC" = "#F39B7F"
)

# Create dot plot with mean bar - all glycoproteins by location (6 categories)
Figure2D <- ggplot(Figure2D_plot_data, aes(x = Location_plot, y = Ratio_144_138, color = Location_plot)) +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.4, color = "black") +
  geom_text(
    data = sample_sizes,
    aes(x = Location_plot, y = 0.02, label = paste0("n=", n)),
    size = 1.8, vjust = 0, color = "black"
  ) +
  scale_color_manual(values = colors_location) +
  coord_cartesian(ylim = c(0, 0.9), clip = "off") +
  scale_y_continuous(breaks = seq(0, 0.9, 0.3)) +
  labs(
    x = NULL,
    y = "Ratio 144/138"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica", color = "black"),
    axis.text = element_text(color = "black", size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 8),
    legend.position = "none",
    plot.margin = margin(10, 5, 5, 5)
  )

# Save plot
ggsave(
  filename = paste0(figure_file_path, "Figure2/Figure2D.pdf"),
  plot = Figure2D,
  width = 2, height = 2, units = "in"
)

cat("\nFigure 2D saved to:", figure_file_path, "Figure2/Figure2D.pdf\n")

