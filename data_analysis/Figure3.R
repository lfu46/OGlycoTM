# Figure 3: O-GlcNAc Differential Analysis Visualization

library(tidyverse)
library(ggpubr)
library(rstatix)

# Source data paths, colors, and differential analysis results
source('data_source_DE.R')

# Figure 3A ---------------------------------------------------------------
# Violin boxplot comparing O-GlcNAc protein logFC distribution across cell types
# Statistical test: Kolmogorov-Smirnov test (compares distribution shapes)

# Combine logFC results from all cell types
OGlcNAc_protein_logFC_HEK293T <- OGlcNAc_protein_DE_HEK293T %>%
  mutate(CellType = "HEK293T") %>%
  select(Protein.ID, logFC, CellType)

OGlcNAc_protein_logFC_HepG2 <- OGlcNAc_protein_DE_HepG2 %>%
  mutate(CellType = "HepG2") %>%
  select(Protein.ID, logFC, CellType)

OGlcNAc_protein_logFC_Jurkat <- OGlcNAc_protein_DE_Jurkat %>%
  mutate(CellType = "Jurkat") %>%
  select(Protein.ID, logFC, CellType)

OGlcNAc_protein_logFC_combined <- bind_rows(
  OGlcNAc_protein_logFC_HEK293T,
  OGlcNAc_protein_logFC_HepG2,
  OGlcNAc_protein_logFC_Jurkat
)

# Set factor levels for consistent ordering
OGlcNAc_protein_logFC_combined <- OGlcNAc_protein_logFC_combined %>%
  mutate(CellType = factor(CellType, levels = c("HEK293T", "HepG2", "Jurkat")))

# Kolmogorov-Smirnov tests between cell types
ks_HEK293T_HepG2 <- ks.test(
  OGlcNAc_protein_DE_HEK293T$logFC,
  OGlcNAc_protein_DE_HepG2$logFC
)
ks_HEK293T_Jurkat <- ks.test(
  OGlcNAc_protein_DE_HEK293T$logFC,
  OGlcNAc_protein_DE_Jurkat$logFC
)
ks_HepG2_Jurkat <- ks.test(
  OGlcNAc_protein_DE_HepG2$logFC,
  OGlcNAc_protein_DE_Jurkat$logFC
)

# Print KS test results
cat("Figure 3A - KS Test Results:\n")
cat("HEK293T vs HepG2: p =", format(ks_HEK293T_HepG2$p.value, digits = 4), "\n")
cat("HEK293T vs Jurkat: p =", format(ks_HEK293T_Jurkat$p.value, digits = 4), "\n")
cat("HepG2 vs Jurkat: p =", format(ks_HepG2_Jurkat$p.value, digits = 4), "\n")

# Wilcoxon test using rstatix (for comparison)
OGlcNAc_protein_wilcox_test <- OGlcNAc_protein_logFC_combined %>%
  wilcox_test(logFC ~ CellType) %>%
  add_significance("p")

cat("\nFigure 3A - Wilcoxon Test Results:\n")
print(OGlcNAc_protein_wilcox_test)

# Create significance annotation dataframe (using KS test)
OGlcNAc_protein_ks_test <- tribble(
  ~.y., ~group1, ~group2, ~p,
  "logFC", "HEK293T", "HepG2", ks_HEK293T_HepG2$p.value,
  "logFC", "HEK293T", "Jurkat", ks_HEK293T_Jurkat$p.value,
  "logFC", "HepG2", "Jurkat", ks_HepG2_Jurkat$p.value
) %>%
  add_significance("p")

# Create violin boxplot
Figure3A <- OGlcNAc_protein_logFC_combined %>%
  ggplot(aes(x = CellType, y = logFC)) +
  geom_violin(aes(fill = CellType), color = "transparent") +
  geom_boxplot(color = "black", outliers = FALSE, width = 0.2, linewidth = 0.3) +
  scale_fill_manual(values = colors_cell) +
  labs(
    x = "",
    y = expression(log[2]*"(Tuni/Ctrl)")
  ) +
  stat_pvalue_manual(
    data = OGlcNAc_protein_ks_test,
    label = "p.signif",
    tip.length = 0,
    size = 3,
    y.position = c(1.5, 1.9, 1.7)
  ) +
  coord_cartesian(ylim = c(-2, 2.3)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 9, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.position = "none"
  )

print(Figure3A)

ggsave(
  filename = paste0(figure_file_path, "Figure3/Figure3A.pdf"),
  plot = Figure3A,
  width = 1.2, height = 1.5, units = "in"
)

cat("\nFigure 3A saved to:", figure_file_path, "Figure3/\n")

# Figure 3B ---------------------------------------------------------------
# Split violin plot comparing O-GlcNAc protein vs whole proteome abundance changes
# Statistical test: Kolmogorov-Smirnov test

library(introdataviz)

# Extract information for HEK293T
OGlcNAc_WP_distribution_HEK293T <- OGlcNAc_protein_DE_HEK293T |>
  select(Protein.ID, logFC_OGlcNAc = logFC) |>
  left_join(WP_protein_DE_HEK293T, by = join_by(Protein.ID == UniProt_Accession)) |>
  select(Protein.ID, logFC_OGlcNAc, logFC_WP = logFC) |>
  filter(!is.na(logFC_WP)) |>
  pivot_longer(cols = logFC_OGlcNAc:logFC_WP, names_to = "Exp", values_to = "logFC") |>
  mutate(Cell = "HEK293T")

# Extract information for HepG2
OGlcNAc_WP_distribution_HepG2 <- OGlcNAc_protein_DE_HepG2 |>
  select(Protein.ID, logFC_OGlcNAc = logFC) |>
  left_join(WP_protein_DE_HepG2, by = join_by(Protein.ID == UniProt_Accession)) |>
  select(Protein.ID, logFC_OGlcNAc, logFC_WP = logFC) |>
  filter(!is.na(logFC_WP)) |>
  pivot_longer(cols = logFC_OGlcNAc:logFC_WP, names_to = "Exp", values_to = "logFC") |>
  mutate(Cell = "HepG2")

# Extract information for Jurkat
OGlcNAc_WP_distribution_Jurkat <- OGlcNAc_protein_DE_Jurkat |>
  select(Protein.ID, logFC_OGlcNAc = logFC) |>
  left_join(WP_protein_DE_Jurkat, by = join_by(Protein.ID == UniProt_Accession)) |>
  select(Protein.ID, logFC_OGlcNAc, logFC_WP = logFC) |>
  filter(!is.na(logFC_WP)) |>
  pivot_longer(cols = logFC_OGlcNAc:logFC_WP, names_to = "Exp", values_to = "logFC") |>
  mutate(Cell = "Jurkat")

# Combine results for all cell types
OGlcNAc_WP_distribution_combined <- bind_rows(
  OGlcNAc_WP_distribution_HEK293T,
  OGlcNAc_WP_distribution_HepG2,
  OGlcNAc_WP_distribution_Jurkat
)

# Print sample sizes
cat("\nFigure 3B - Sample sizes (matched proteins):\n")
cat("HEK293T:", nrow(OGlcNAc_WP_distribution_HEK293T) / 2, "\n")
cat("HepG2:", nrow(OGlcNAc_WP_distribution_HepG2) / 2, "\n")
cat("Jurkat:", nrow(OGlcNAc_WP_distribution_Jurkat) / 2, "\n")

# KS tests
ks_HEK293T <- ks.test(
  OGlcNAc_WP_distribution_HEK293T |> filter(Exp == "logFC_OGlcNAc") |> pull(logFC),
  OGlcNAc_WP_distribution_HEK293T |> filter(Exp == "logFC_WP") |> pull(logFC)
)

ks_HepG2 <- ks.test(
  OGlcNAc_WP_distribution_HepG2 |> filter(Exp == "logFC_OGlcNAc") |> pull(logFC),
  OGlcNAc_WP_distribution_HepG2 |> filter(Exp == "logFC_WP") |> pull(logFC)
)

ks_Jurkat <- ks.test(
  OGlcNAc_WP_distribution_Jurkat |> filter(Exp == "logFC_OGlcNAc") |> pull(logFC),
  OGlcNAc_WP_distribution_Jurkat |> filter(Exp == "logFC_WP") |> pull(logFC)
)

cat("\nFigure 3B - KS Test Results (O-GlcNAc vs WP):\n")
cat("HEK293T: p =", format(ks_HEK293T$p.value, digits = 4), "\n")
cat("HepG2: p =", format(ks_HepG2$p.value, digits = 4), "\n")
cat("Jurkat: p =", format(ks_Jurkat$p.value, digits = 4), "\n")

# Function to convert p-value to significance label
get_signif_label <- function(p) {
  if (p < 0.0001) return("****")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

# Create KS test label dataframe for annotation
ks_test_label <- tibble(
  Cell = c("HEK293T", "HepG2", "Jurkat"),
  p_signif = c(
    get_signif_label(ks_HEK293T$p.value),
    get_signif_label(ks_HepG2$p.value),
    get_signif_label(ks_Jurkat$p.value)
  ),
  logFC = 1.5
) |>
  mutate(Cell = factor(Cell, levels = c("HEK293T", "HepG2", "Jurkat")))

# Prepare data for split violin plot (without faceting)
Figure3B_data <- OGlcNAc_WP_distribution_combined |>
  mutate(
    Cell = factor(Cell, levels = c("HEK293T", "HepG2", "Jurkat")),
    # Create unique fill variable combining cell type and experiment
    fill_group = case_when(
      Exp == "logFC_OGlcNAc" & Cell == "HEK293T" ~ "OGlcNAc_HEK293T",
      Exp == "logFC_OGlcNAc" & Cell == "HepG2" ~ "OGlcNAc_HepG2",
      Exp == "logFC_OGlcNAc" & Cell == "Jurkat" ~ "OGlcNAc_Jurkat",
      Exp == "logFC_WP" ~ "WP"
    ),
    fill_group = factor(fill_group, levels = c("OGlcNAc_HEK293T", "OGlcNAc_HepG2", "OGlcNAc_Jurkat", "WP"))
  )

# Split violin plot (combined, no faceting)
Figure3B <- Figure3B_data |>
  ggplot(aes(x = Cell, y = logFC, fill = fill_group)) +
  geom_split_violin(color = "transparent") +
  geom_text(data = ks_test_label, aes(x = Cell, y = logFC, label = p_signif),
            inherit.aes = FALSE, size = 3) +
  scale_fill_manual(values = c(
    "OGlcNAc_HEK293T" = unname(colors_cell["HEK293T"]),
    "OGlcNAc_HepG2" = unname(colors_cell["HepG2"]),
    "OGlcNAc_Jurkat" = unname(colors_cell["Jurkat"]),
    "WP" = "gray70"
  )) +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)"), fill = "") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 9, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.position = "none"
  )

print(Figure3B)

ggsave(
  filename = paste0(figure_file_path, "Figure3/Figure3B.pdf"),
  plot = Figure3B,
  width = 1.2, height = 1.5, units = "in"
)

cat("\nFigure 3B saved to:", figure_file_path, "Figure3/\n")

# Figure 3C ---------------------------------------------------------------
# UMAP plot of overlapping O-GlcNAc proteins using logFC as features
# Highlights commonly up/downregulated proteins across cell types

library(reticulate)

# Configure Python environment
options(reticulate.conda_binary = "/opt/anaconda3/bin/conda")
use_python("/opt/anaconda3/envs/umap_env/bin/python", required = TRUE)

# Source Python UMAP script
source_python("umap_analysis.py")

# Find overlapping proteins across all 3 cell types
overlapping_proteins <- Reduce(intersect, list(
  OGlcNAc_protein_DE_HEK293T$Protein.ID,
  OGlcNAc_protein_DE_HepG2$Protein.ID,
  OGlcNAc_protein_DE_Jurkat$Protein.ID
))

cat("Figure 3C - Number of overlapping proteins:", length(overlapping_proteins), "\n")

# Create logFC matrix for UMAP (proteins as rows, cell types as columns)
logFC_matrix <- tibble(Protein.ID = overlapping_proteins) %>%
  left_join(
    OGlcNAc_protein_DE_HEK293T %>% select(Protein.ID, logFC) %>% rename(logFC_HEK293T = logFC),
    by = "Protein.ID"
  ) %>%
  left_join(
    OGlcNAc_protein_DE_HepG2 %>% select(Protein.ID, logFC) %>% rename(logFC_HepG2 = logFC),
    by = "Protein.ID"
  ) %>%
  left_join(
    OGlcNAc_protein_DE_Jurkat %>% select(Protein.ID, logFC) %>% rename(logFC_Jurkat = logFC),
    by = "Protein.ID"
  )

# Prepare data for Python UMAP
protein_ids <- logFC_matrix$Protein.ID
logFC_data <- logFC_matrix %>% select(starts_with("logFC_"))

# Run UMAP via Python
umap_result <- run_umap_proteins(
  logfc_matrix = logFC_data,
  protein_ids = protein_ids,
  n_neighbors = 15L,
  min_dist = 0.1,
  random_state = 42L
)

# Convert to tibble
umap_df <- as_tibble(umap_result)

# Classify proteins as commonly up, commonly down, or other
umap_df <- umap_df %>%
  mutate(
    Regulation = case_when(
      Protein.ID %in% OGlcNAc_protein_commonly_up$Protein.ID ~ "Commonly Up",
      Protein.ID %in% OGlcNAc_protein_commonly_down$Protein.ID ~ "Commonly Down",
      TRUE ~ "Other"
    )
  ) %>%
  mutate(Regulation = factor(Regulation, levels = c("Commonly Up", "Commonly Down", "Other")))

# Print counts
cat("Protein classification:\n")
print(table(umap_df$Regulation))

# Define colors and shapes
colors_regulation <- c(
  "Commonly Up" = "#F39B7F",
  "Commonly Down" = "#4DBBD5",
  "Other" = "grey70"
)

shapes_regulation <- c(
  "Commonly Up" = 17,    # solid triangle up
  "Commonly Down" = 15,  # solid square
  "Other" = 16           # solid circle
)

# Create UMAP plot (plot "Other" first, then highlighted proteins on top)
Figure3C <- ggplot() +
  # Plot "Other" proteins first (background)
  geom_point(
    data = umap_df %>% filter(Regulation == "Other"),
    aes(x = UMAP1, y = UMAP2, color = Regulation, shape = Regulation),
    size = 1, alpha = 0.5
  ) +
  # Plot commonly regulated proteins on top
  geom_point(
    data = umap_df %>% filter(Regulation != "Other"),
    aes(x = UMAP1, y = UMAP2, color = Regulation, shape = Regulation),
    size = 2, alpha = 0.8
  ) +
  scale_color_manual(
    values = colors_regulation,
    breaks = c("Commonly Up", "Commonly Down")
  ) +
  scale_shape_manual(
    values = shapes_regulation,
    breaks = c("Commonly Up", "Commonly Down")
  ) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    color = "",
    shape = ""
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 9, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.2, "cm")
  )

print(Figure3C)

ggsave(
  filename = paste0(figure_file_path, "Figure3/Figure3C.pdf"),
  plot = Figure3C,
  width = 1.5, height = 1.8, units = "in"
)

cat("\nFigure 3C saved to:", figure_file_path, "Figure3/\n")

# Figure 3D ---------------------------------------------------------------
# GO enrichment dot plot for commonly upregulated O-GlcNAc proteins

# library(clusterProfiler)
# library(org.Hs.eg.db)

# # Get overlapping proteins (if not already created in Figure 3C)
# if (!exists("overlapping_proteins")) {
#   overlapping_proteins <- Reduce(intersect, list(
#     OGlcNAc_protein_DE_HEK293T$Protein.ID,
#     OGlcNAc_protein_DE_HepG2$Protein.ID,
#     OGlcNAc_protein_DE_Jurkat$Protein.ID
#   ))
# }

# # Get protein lists for enrichment
# commonly_up_proteins <- OGlcNAc_protein_commonly_up$Protein.ID
# commonly_down_proteins <- OGlcNAc_protein_commonly_down$Protein.ID

# cat("Figure 3D - Commonly upregulated proteins:", length(commonly_up_proteins), "\n")
# cat("Figure 3D - Commonly downregulated proteins:", length(commonly_down_proteins), "\n")

# # Background 1: Overlapping O-GlcNAc proteins (402)
# background_overlapping <- overlapping_proteins
# cat("Background (overlapping):", length(background_overlapping), "proteins\n")

# # Background 2: Total quantified O-GlcNAc proteins (union of all 3 cell types)
# background_total <- Reduce(union, list(
#   OGlcNAc_protein_DE_HEK293T$Protein.ID,
#   OGlcNAc_protein_DE_HepG2$Protein.ID,
#   OGlcNAc_protein_DE_Jurkat$Protein.ID
# ))
# cat("Background (total quantified):", length(background_total), "proteins\n")

# # --- Enrichment with overlapping background ---

# # GO enrichment for commonly upregulated proteins
# commonly_up_GO_overlapping <- enrichGO(
#   gene = commonly_up_proteins,
#   OrgDb = org.Hs.eg.db,
#   universe = background_overlapping,
#   keyType = 'UNIPROT',
#   ont = 'ALL',
#   pvalueCutoff = 1,
#   qvalueCutoff = 1
# )

# write_csv(
#   commonly_up_GO_overlapping@result,
#   paste0(source_file_path, 'enrichment/Figure3D_commonly_up_GO_overlapping.csv')
# )

# # GO enrichment for commonly downregulated proteins
# commonly_down_GO_overlapping <- enrichGO(
#   gene = commonly_down_proteins,
#   OrgDb = org.Hs.eg.db,
#   universe = background_overlapping,
#   keyType = 'UNIPROT',
#   ont = 'ALL',
#   pvalueCutoff = 1,
#   qvalueCutoff = 1
# )

# write_csv(
#   commonly_down_GO_overlapping@result,
#   paste0(source_file_path, 'enrichment/Figure3D_commonly_down_GO_overlapping.csv')
# )

# # --- Enrichment with total O-GlcNAc background ---

# # GO enrichment for commonly upregulated proteins
# commonly_up_GO_total <- enrichGO(
#   gene = commonly_up_proteins,
#   OrgDb = org.Hs.eg.db,
#   universe = background_total,
#   keyType = 'UNIPROT',
#   ont = 'ALL',
#   pvalueCutoff = 1,
#   qvalueCutoff = 1
# )

# write_csv(
#   commonly_up_GO_total@result,
#   paste0(source_file_path, 'enrichment/Figure3D_commonly_up_GO_total.csv')
# )

# # GO enrichment for commonly downregulated proteins
# commonly_down_GO_total <- enrichGO(
#   gene = commonly_down_proteins,
#   OrgDb = org.Hs.eg.db,
#   universe = background_total,
#   keyType = 'UNIPROT',
#   ont = 'ALL',
#   pvalueCutoff = 1,
#   qvalueCutoff = 1
# )

# write_csv(
#   commonly_down_GO_total@result,
#   paste0(source_file_path, 'enrichment/Figure3D_commonly_down_GO_total.csv')
# )

# cat("\nEnrichment results saved to:", source_file_path, "enrichment/\n")

# Select GO terms for commonly upregulated proteins (manually select 3-4 terms)
commonly_up_GO_selected <- read_csv(
  paste0(source_file_path, 'enrichment/Figure3E_commonly_up_GO_total.csv'),
  show_col_types = FALSE
) |>
  filter(
    Description %in% c(
      'organophosphate metabolic process',
      'regulation of translation',
      'response to glucose',
      'positive regulation of RNA splicing'
    )
  ) |>
  dplyr::select(Description, ONTOLOGY, pvalue, Count) |>
  mutate(log_pvalue = -log10(pvalue))

# Select GO terms for commonly downregulated proteins (manually select 3-4 terms)
commonly_down_GO_selected <- read_csv(
  paste0(source_file_path, 'enrichment/Figure3E_commonly_down_GO_total.csv'),
  show_col_types = FALSE
) |>
  filter(
    Description %in% c(
      'nuclear protein-containing complex',
      'regulation of response to stress',
      'regulation of translational initiation',
      'cytoplasmic stress granule'
    )
  ) |>
  dplyr::select(Description, ONTOLOGY, pvalue, Count) |>
  mutate(log_pvalue = -log10(pvalue))

# Create dot plot for commonly upregulated proteins (no legend)
Figure3D_up <- ggplot(commonly_up_GO_selected, aes(x = log_pvalue, y = reorder(Description, log_pvalue))) +
  geom_segment(
    aes(x = 0, xend = log_pvalue, y = Description, yend = Description),
    linetype = "dashed",
    color = "grey70",
    linewidth = 0.3
  ) +
  geom_point(aes(size = Count), color = "#F39B7F") +
  scale_size_continuous(range = c(2, 5)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = expression(-log[10]('p value')),
    y = NULL,
    size = "Count"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 9),
    axis.text.y = element_text(lineheight = 0.8),
    axis.title = element_text(size = 9),
    legend.position = "none"
  )

# Create dot plot for commonly downregulated proteins (with legend)
Figure3D_down <- ggplot(commonly_down_GO_selected, aes(x = log_pvalue, y = reorder(Description, log_pvalue))) +
  geom_segment(
    aes(x = 0, xend = log_pvalue, y = Description, yend = Description),
    linetype = "dashed",
    color = "grey70",
    linewidth = 0.3
  ) +
  geom_point(aes(size = Count), color = "#4DBBD5") +
  scale_size_continuous(range = c(2, 5)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = expression(-log[10]('p value')),
    y = NULL,
    size = "Count"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 9),
    axis.text.y = element_text(lineheight = 0.8),
    axis.title = element_text(size = 9),
    legend.position = "right",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.key.size = unit(0.4, "cm")
  )

# Combine plots horizontally with shared legend
library(patchwork)

Figure3D <- Figure3D_up + Figure3D_down +
  plot_layout(ncol = 2, guides = "collect")

print(Figure3D)

ggsave(
  filename = paste0(figure_file_path, "Figure3/Figure3D.pdf"),
  plot = Figure3D,
  width = 6, height = 1.5, units = "in"
)

cat("\nFigure 3D saved to:", figure_file_path, "Figure3/\n")

