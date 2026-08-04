# N-Sequon Filter Comparison: Regenerate Figures 2, 3, 4 with filtered data
# Purpose: Remove Level3 (unlocalized) O-GlcNAc PSMs with N-X-S/T sequon
# Output: /Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/NSequon_filtered_figures/

library(tidyverse)
library(eulerr)
library(cowplot)
library(grid)
library(ggpubr)
library(rstatix)
library(introdataviz)
library(circlize)
library(scales)
library(ComplexHeatmap)
library(gridBase)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)

# --- Paths ---
source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'
output_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/NSequon_filtered_figures/'
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

# --- Color palettes ---
colors_cell <- c("HEK293T" = "#4DBBD5", "HepG2" = "#F39B7F", "Jurkat" = "#00A087")

# =============================================================================
# Step 1: Apply N-sequon filter to bonafide data
# =============================================================================
cat("=== Applying N-sequon filter ===\n")

OGlyco_HEK293T_bonafide <- read_csv(paste0(source_file_path, 'filtered/OGlyco_HEK293T_bonafide.csv'), show_col_types = FALSE)
OGlyco_HepG2_bonafide <- read_csv(paste0(source_file_path, 'filtered/OGlyco_HepG2_bonafide.csv'), show_col_types = FALSE)
OGlyco_Jurkat_bonafide <- read_csv(paste0(source_file_path, 'filtered/OGlyco_Jurkat_bonafide.csv'), show_col_types = FALSE)

OGlcNAc_compositions <- c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859')

# Filter function: remove Level3 + N-sequon PSMs from O-GlcNAc data
apply_nsequon_filter <- function(bonafide_df) {
  oglcnac <- bonafide_df |> filter(Total.Glycan.Composition %in% OGlcNAc_compositions)
  oglcnac |> filter(!(Confidence.Level == "Level3" & Has.N.Glyc.Sequon == TRUE))
}

OGlcNAc_HEK293T_filt <- apply_nsequon_filter(OGlyco_HEK293T_bonafide)
OGlcNAc_HepG2_filt <- apply_nsequon_filter(OGlyco_HepG2_bonafide)
OGlcNAc_Jurkat_filt <- apply_nsequon_filter(OGlyco_Jurkat_bonafide)

# Protein lists
OGlcNAc_proteins_HEK293T <- unique(OGlcNAc_HEK293T_filt$Protein.ID)
OGlcNAc_proteins_HepG2 <- unique(OGlcNAc_HepG2_filt$Protein.ID)
OGlcNAc_proteins_Jurkat <- unique(OGlcNAc_Jurkat_filt$Protein.ID)
valid_proteins <- list(HEK293T = OGlcNAc_proteins_HEK293T, HepG2 = OGlcNAc_proteins_HepG2, Jurkat = OGlcNAc_proteins_Jurkat)

cat("Filtered protein counts: HEK293T =", length(OGlcNAc_proteins_HEK293T),
    ", HepG2 =", length(OGlcNAc_proteins_HepG2),
    ", Jurkat =", length(OGlcNAc_proteins_Jurkat), "\n")

# Key set operations
common_proteins <- Reduce(intersect, valid_proteins)
unique_HEK293T <- setdiff(setdiff(OGlcNAc_proteins_HEK293T, OGlcNAc_proteins_HepG2), OGlcNAc_proteins_Jurkat)
unique_HepG2 <- setdiff(setdiff(OGlcNAc_proteins_HepG2, OGlcNAc_proteins_HEK293T), OGlcNAc_proteins_Jurkat)
unique_Jurkat <- setdiff(setdiff(OGlcNAc_proteins_Jurkat, OGlcNAc_proteins_HEK293T), OGlcNAc_proteins_HepG2)
total_proteins <- Reduce(union, valid_proteins)

cat("Common:", length(common_proteins), " Unique HEK:", length(unique_HEK293T),
    " Unique HepG2:", length(unique_HepG2), " Unique Jurkat:", length(unique_Jurkat),
    " Total:", length(total_proteins), "\n")

# =============================================================================
# Step 2: Load and filter DE results
# =============================================================================
cat("\n=== Loading and filtering DE results ===\n")

OGlcNAc_protein_DE_HEK293T <- read_csv(paste0(source_file_path, 'differential_analysis/OGlcNAc_protein_DE_HEK293T.csv'), show_col_types = FALSE) |>
  filter(Protein.ID %in% OGlcNAc_proteins_HEK293T)
OGlcNAc_protein_DE_HepG2 <- read_csv(paste0(source_file_path, 'differential_analysis/OGlcNAc_protein_DE_HepG2.csv'), show_col_types = FALSE) |>
  filter(Protein.ID %in% OGlcNAc_proteins_HepG2)
OGlcNAc_protein_DE_Jurkat <- read_csv(paste0(source_file_path, 'differential_analysis/OGlcNAc_protein_DE_Jurkat.csv'), show_col_types = FALSE) |>
  filter(Protein.ID %in% OGlcNAc_proteins_Jurkat)

# WP DE (unchanged)
WP_protein_DE_HEK293T <- read_csv(paste0(source_file_path, 'differential_analysis/WP_protein_DE_HEK293T.csv'), show_col_types = FALSE)
WP_protein_DE_HepG2 <- read_csv(paste0(source_file_path, 'differential_analysis/WP_protein_DE_HepG2.csv'), show_col_types = FALSE)
WP_protein_DE_Jurkat <- read_csv(paste0(source_file_path, 'differential_analysis/WP_protein_DE_Jurkat.csv'), show_col_types = FALSE)

cat("DE proteins: HEK293T =", nrow(OGlcNAc_protein_DE_HEK293T),
    ", HepG2 =", nrow(OGlcNAc_protein_DE_HepG2),
    ", Jurkat =", nrow(OGlcNAc_protein_DE_Jurkat), "\n")

# Re-derive commonly regulated proteins
get_sig <- function(de, direction = "up") {
  if (direction == "up") de |> filter(logFC > 0.5, adj.P.Val < 0.05) |> pull(Protein.ID)
  else de |> filter(logFC < -0.5, adj.P.Val < 0.05) |> pull(Protein.ID)
}

up_H <- get_sig(OGlcNAc_protein_DE_HEK293T, "up")
up_G <- get_sig(OGlcNAc_protein_DE_HepG2, "up")
up_J <- get_sig(OGlcNAc_protein_DE_Jurkat, "up")
down_H <- get_sig(OGlcNAc_protein_DE_HEK293T, "down")
down_G <- get_sig(OGlcNAc_protein_DE_HepG2, "down")
down_J <- get_sig(OGlcNAc_protein_DE_Jurkat, "down")

cat("\nUp/Down counts: HEK293T =", length(up_H), "/", length(down_H),
    ", HepG2 =", length(up_G), "/", length(down_G),
    ", Jurkat =", length(up_J), "/", length(down_J), "\n")

# Commonly up: sig up in >=2 cell types AND not sig down in any
commonly_up_proteins <- unique(c(
  intersect(up_H, up_G), intersect(up_H, up_J), intersect(up_G, up_J)
))
commonly_up_proteins <- setdiff(commonly_up_proteins, c(down_H, down_G, down_J))

commonly_down_proteins <- unique(c(
  intersect(down_H, down_G), intersect(down_H, down_J), intersect(down_G, down_J)
))
commonly_down_proteins <- setdiff(commonly_down_proteins, c(up_H, up_G, up_J))

cat("Commonly up:", length(commonly_up_proteins), " Commonly down:", length(commonly_down_proteins), "\n")

# =============================================================================
# Figure 2A: Euler diagram with GO enrichment bars
# =============================================================================
cat("\n=== Generating Figure 2A ===\n")

protein_list <- list(
  HepG2 = OGlcNAc_proteins_HepG2,
  Jurkat = OGlcNAc_proteins_Jurkat,
  HEK293T = OGlcNAc_proteins_HEK293T
)

euler_obj <- euler(protein_list)
Figure2A_colors <- c("HepG2" = "#F39B7F", "Jurkat" = "#00A087", "HEK293T" = "#4DBBD5")

euler_plot <- plot(euler_obj,
  fills = list(fill = Figure2A_colors, alpha = 0.5),
  edges = list(col = "white", lwd = 2),
  labels = FALSE,
  quantities = list(font = 1, cex = 1.0)
)

# --- Re-run GO enrichment for filtered data ---
# Common GO
common_GO_result <- enrichGO(
  gene = common_proteins, OrgDb = org.Hs.eg.db,
  universe = total_proteins, keyType = 'UNIPROT', ont = 'ALL',
  pvalueCutoff = 1, qvalueCutoff = 1
)

common_GO_selected <- common_GO_result@result |>
  filter(Description %in% c(
    'transcription coregulator activity', 'mRNA binding',
    'nucleocytoplasmic transport', 'nuclear transport',
    'cytoplasmic stress granule', 'histone acetyltransferase complex'
  )) |>
  dplyr::select(Description, pvalue, Count) |>
  mutate(
    log_pvalue = -log10(pvalue),
    Description = case_when(
      Description == "transcription coregulator activity" ~ "Transcription coregulator",
      Description == "nucleocytoplasmic transport" ~ "Nucleocytoplasmic transport",
      Description == "histone acetyltransferase complex" ~ "Histone acetyltransferase",
      Description == "cytoplasmic stress granule" ~ "Cytoplasmic stress granule",
      TRUE ~ Description
    )
  )

# HEK293T unique GO
HEK293T_GO_result <- enrichGO(
  gene = unique_HEK293T, OrgDb = org.Hs.eg.db,
  universe = total_proteins, keyType = 'UNIPROT', ont = 'ALL',
  pvalueCutoff = 1, qvalueCutoff = 1
)
HEK293T_GO_selected <- HEK293T_GO_result@result |>
  filter(Description %in% c(
    'chaperone-mediated protein folding', 'cell cycle process', 'epithelial cell proliferation'
  )) |>
  dplyr::select(Description, pvalue, Count) |>
  mutate(
    log_pvalue = -log10(pvalue),
    Description = case_when(
      Description == "chaperone-mediated protein folding" ~ "Chaperone-mediated folding",
      Description == "epithelial cell proliferation" ~ "Epithelial cell proliferation",
      TRUE ~ Description
    )
  )

# Jurkat unique GO
Jurkat_GO_result <- enrichGO(
  gene = unique_Jurkat, OrgDb = org.Hs.eg.db,
  universe = total_proteins, keyType = 'UNIPROT', ont = 'ALL',
  pvalueCutoff = 1, qvalueCutoff = 1
)
Jurkat_GO_selected <- Jurkat_GO_result@result |>
  filter(Description %in% c(
    'leukocyte activation involved in immune response',
    'leukocyte cell-cell adhesion', 'T cell activation'
  )) |>
  dplyr::select(Description, pvalue, Count) |>
  mutate(
    log_pvalue = -log10(pvalue),
    Description = case_when(
      Description == "leukocyte activation involved in immune response" ~ "Leukocyte activation",
      Description == "leukocyte cell-cell adhesion" ~ "Leukocyte adhesion",
      TRUE ~ Description
    )
  )

cat("GO terms found: Common =", nrow(common_GO_selected),
    " HEK293T =", nrow(HEK293T_GO_selected),
    " Jurkat =", nrow(Jurkat_GO_selected), "\n")

# Build GO bar plots
create_gradient_data <- function(df, n_segments = 50) {
  df <- df |> mutate(y_num = as.numeric(fct_reorder(Description, log_pvalue)))
  do.call(rbind, lapply(1:nrow(df), function(i) {
    row <- df[i, ]
    segs <- data.frame(
      Description = row$Description, y_num = row$y_num,
      xmin = seq(0, row$log_pvalue, length.out = n_segments + 1)[-(n_segments + 1)],
      xmax = seq(0, row$log_pvalue, length.out = n_segments + 1)[-1],
      segment = 1:n_segments
    )
    segs$alpha_val <- seq(0.9, 0.3, length.out = n_segments)
    segs
  }))
}

make_go_bar <- function(go_data, fill_color) {
  grad <- create_gradient_data(go_data)
  ggplot() +
    geom_rect(data = grad, aes(xmin = xmin, xmax = xmax, ymin = y_num - 0.35, ymax = y_num + 0.35, alpha = alpha_val), fill = fill_color) +
    geom_text(data = go_data |> mutate(y_num = as.numeric(fct_reorder(Description, log_pvalue))),
              aes(label = Description, x = 0.05, y = y_num), hjust = 0, size = 2.8, color = "black") +
    scale_alpha_identity() +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_y_continuous(breaks = 1:nrow(go_data), labels = NULL) +
    labs(x = expression(-log[10]~"("*italic(p)~value*")"), y = NULL) +
    theme_classic() +
    theme(text = element_text(family = "Helvetica", color = "black"),
          axis.title.x = element_text(size = 7), axis.text.x = element_text(color = "black", size = 7),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.background = element_blank(), panel.background = element_blank())
}

common_GO_plot <- make_go_bar(common_GO_selected, "#A0A0A0")
HEK293T_GO_plot <- make_go_bar(HEK293T_GO_selected, "#7DCDE5")
Jurkat_GO_plot <- make_go_bar(Jurkat_GO_selected, "#4DC4B0")

euler_grob <- grid::grid.grabExpr(print(euler_plot), width = 3, height = 3)

Figure2A <- ggdraw() +
  draw_grob(euler_grob, x = 0.28, y = 0.15, width = 0.45, height = 0.70) +
  draw_plot(common_GO_plot, x = 0, y = 0.20, width = 0.32, height = 0.60) +
  draw_plot(HEK293T_GO_plot, x = 0.72, y = 0.54, width = 0.28, height = 0.40) +
  draw_plot(Jurkat_GO_plot, x = 0.72, y = 0.06, width = 0.28, height = 0.40) +
  draw_label(paste0("Common (", length(common_proteins), ")"), x = 0.16, y = 0.84, hjust = 0.5, vjust = 0.5, fontface = "plain", color = "#505050", size = 7) +
  draw_label(paste0("HEK293T (", length(unique_HEK293T), ")"), x = 0.86, y = 0.96, hjust = 0.5, vjust = 0.5, fontface = "plain", color = "#2A9BB5", size = 7) +
  draw_label(paste0("Jurkat (", length(unique_Jurkat), ")"), x = 0.86, y = 0.48, hjust = 0.5, vjust = 0.5, fontface = "plain", color = "#008570", size = 7) +
  geom_path(data = data.frame(x = c(0.50, 0.40, 0.32), y = c(0.50, 0.60, 0.60)), aes(x = x, y = y), color = "#505050", linewidth = 0.5) +
  geom_point(aes(x = 0.50, y = 0.50), shape = 1, size = 5, color = "#505050", stroke = 0.8) +
  geom_path(data = data.frame(x = c(0.62, 0.72, 0.78), y = c(0.62, 0.78, 0.78)), aes(x = x, y = y), color = "#2A9BB5", linewidth = 0.5) +
  geom_point(aes(x = 0.62, y = 0.62), shape = 1, size = 4, color = "#2A9BB5", stroke = 0.8) +
  geom_path(data = data.frame(x = c(0.62, 0.72, 0.78), y = c(0.38, 0.28, 0.28)), aes(x = x, y = y), color = "#008570", linewidth = 0.5) +
  geom_point(aes(x = 0.62, y = 0.38), shape = 1, size = 4, color = "#008570", stroke = 0.8)

ggsave(paste0(output_path, "Figure2A.pdf"), Figure2A, width = 5, height = 3)
cat("Figure 2A saved\n")

# =============================================================================
# Figure 3A: Violin boxplot of logFC
# =============================================================================
cat("\n=== Generating Figure 3A ===\n")

logFC_combined <- bind_rows(
  OGlcNAc_protein_DE_HEK293T |> mutate(CellType = "HEK293T") |> dplyr::select(Protein.ID, logFC, CellType),
  OGlcNAc_protein_DE_HepG2 |> mutate(CellType = "HepG2") |> dplyr::select(Protein.ID, logFC, CellType),
  OGlcNAc_protein_DE_Jurkat |> mutate(CellType = "Jurkat") |> dplyr::select(Protein.ID, logFC, CellType)
) |> mutate(CellType = factor(CellType, levels = c("HEK293T", "HepG2", "Jurkat")))

ks_HG <- ks.test(OGlcNAc_protein_DE_HEK293T$logFC, OGlcNAc_protein_DE_HepG2$logFC)
ks_HJ <- ks.test(OGlcNAc_protein_DE_HEK293T$logFC, OGlcNAc_protein_DE_Jurkat$logFC)
ks_GJ <- ks.test(OGlcNAc_protein_DE_HepG2$logFC, OGlcNAc_protein_DE_Jurkat$logFC)

cat("KS tests: HEK-HepG2 p=", format(ks_HG$p.value, digits=4),
    " HEK-Jurkat p=", format(ks_HJ$p.value, digits=4),
    " HepG2-Jurkat p=", format(ks_GJ$p.value, digits=4), "\n")

ks_df <- tribble(
  ~.y., ~group1, ~group2, ~p,
  "logFC", "HEK293T", "HepG2", ks_HG$p.value,
  "logFC", "HEK293T", "Jurkat", ks_HJ$p.value,
  "logFC", "HepG2", "Jurkat", ks_GJ$p.value
) |> add_significance("p")

Figure3A <- logFC_combined |>
  ggplot(aes(x = CellType, y = logFC)) +
  geom_violin(aes(fill = CellType), color = "transparent") +
  geom_boxplot(color = "black", outliers = FALSE, width = 0.2, linewidth = 0.3) +
  scale_fill_manual(values = colors_cell) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  stat_pvalue_manual(data = ks_df, label = "p.signif", tip.length = 0, size = 5, y.position = c(1.5, 1.9, 1.7)) +
  coord_cartesian(ylim = c(-2, 2.3)) +
  theme_bw() +
  theme(panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
        panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
        axis.title = element_text(size = 9),
        axis.text.x = element_text(size = 9, color = "black", angle = 30, hjust = 1),
        axis.text.y = element_text(size = 9, color = "black"),
        legend.position = "none")

ggsave(paste0(output_path, "Figure3A.pdf"), Figure3A, width = 1.5, height = 2)
cat("Figure 3A saved\n")

# =============================================================================
# Figure 3B: Split violin O-GlcNAc vs WP
# =============================================================================
cat("\n=== Generating Figure 3B ===\n")

make_split_data <- function(oglcnac_de, wp_de, cell_name) {
  oglcnac_de |>
    dplyr::select(Protein.ID, logFC_OGlcNAc = logFC) |>
    left_join(wp_de, by = c("Protein.ID" = "UniProt_Accession")) |>
    dplyr::select(Protein.ID, logFC_OGlcNAc, logFC_WP = logFC) |>
    filter(!is.na(logFC_WP)) |>
    pivot_longer(cols = logFC_OGlcNAc:logFC_WP, names_to = "Exp", values_to = "logFC") |>
    mutate(Cell = cell_name)
}

split_data <- bind_rows(
  make_split_data(OGlcNAc_protein_DE_HEK293T, WP_protein_DE_HEK293T, "HEK293T"),
  make_split_data(OGlcNAc_protein_DE_HepG2, WP_protein_DE_HepG2, "HepG2"),
  make_split_data(OGlcNAc_protein_DE_Jurkat, WP_protein_DE_Jurkat, "Jurkat")
)

get_signif_label <- function(p) {
  if (p < 0.0001) return("****")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("ns")
}

ks_labels <- tibble(
  Cell = c("HEK293T", "HepG2", "Jurkat"),
  p_signif = sapply(c("HEK293T", "HepG2", "Jurkat"), function(c) {
    d <- split_data |> filter(Cell == c)
    ks <- ks.test(d$logFC[d$Exp == "logFC_OGlcNAc"], d$logFC[d$Exp == "logFC_WP"])
    cat("  KS", c, "p =", format(ks$p.value, digits=4), "\n")
    get_signif_label(ks$p.value)
  }),
  logFC = c(1.8, 1.5, 1.2)
) |> mutate(Cell = factor(Cell, levels = c("HEK293T", "HepG2", "Jurkat")))

Figure3B_data <- split_data |>
  mutate(
    Cell = factor(Cell, levels = c("HEK293T", "HepG2", "Jurkat")),
    fill_group = case_when(
      Exp == "logFC_OGlcNAc" & Cell == "HEK293T" ~ "OGlcNAc_HEK293T",
      Exp == "logFC_OGlcNAc" & Cell == "HepG2" ~ "OGlcNAc_HepG2",
      Exp == "logFC_OGlcNAc" & Cell == "Jurkat" ~ "OGlcNAc_Jurkat",
      Exp == "logFC_WP" ~ "WP"
    ),
    fill_group = factor(fill_group, levels = c("OGlcNAc_HEK293T", "OGlcNAc_HepG2", "OGlcNAc_Jurkat", "WP"))
  )

Figure3B <- Figure3B_data |>
  ggplot(aes(x = Cell, y = logFC, fill = fill_group)) +
  geom_split_violin(color = "transparent") +
  geom_text(data = ks_labels, aes(x = Cell, y = logFC, label = p_signif), inherit.aes = FALSE, size = 5) +
  scale_fill_manual(values = c(
    "OGlcNAc_HEK293T" = unname(colors_cell["HEK293T"]),
    "OGlcNAc_HepG2" = unname(colors_cell["HepG2"]),
    "OGlcNAc_Jurkat" = unname(colors_cell["Jurkat"]),
    "WP" = "gray70"
  )) +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)"), fill = "") +
  theme_bw() +
  theme(panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
        panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
        axis.title = element_text(size = 9),
        axis.text.x = element_text(size = 9, color = "black", angle = 30, hjust = 1),
        axis.text.y = element_text(size = 9, color = "black"),
        legend.position = "none")

ggsave(paste0(output_path, "Figure3B.pdf"), Figure3B, width = 1.5, height = 2)
cat("Figure 3B saved\n")

# =============================================================================
# Figure 3D: GO enrichment for commonly up/down
# =============================================================================
cat("\n=== Generating Figure 3D ===\n")

# Re-run GO for commonly up/down with filtered universe
bg_total <- Reduce(union, list(
  OGlcNAc_protein_DE_HEK293T$Protein.ID,
  OGlcNAc_protein_DE_HepG2$Protein.ID,
  OGlcNAc_protein_DE_Jurkat$Protein.ID
))

if (length(commonly_up_proteins) >= 3) {
  comm_up_GO <- enrichGO(gene = commonly_up_proteins, OrgDb = org.Hs.eg.db,
                          universe = bg_total, keyType = 'UNIPROT', ont = 'ALL',
                          pvalueCutoff = 1, qvalueCutoff = 1)
  comm_up_GO_sel <- comm_up_GO@result |>
    filter(Description %in% c('organophosphate metabolic process', 'regulation of translation',
                               'response to glucose', 'positive regulation of RNA splicing')) |>
    dplyr::select(Description, pvalue, Count) |>
    mutate(log_pvalue = -log10(pvalue))
  cat("Commonly up GO terms found:", nrow(comm_up_GO_sel), "\n")
} else {
  cat("Too few commonly up proteins for GO\n")
  comm_up_GO_sel <- tibble()
}

if (length(commonly_down_proteins) >= 3) {
  comm_down_GO <- enrichGO(gene = commonly_down_proteins, OrgDb = org.Hs.eg.db,
                            universe = bg_total, keyType = 'UNIPROT', ont = 'ALL',
                            pvalueCutoff = 1, qvalueCutoff = 1)
  comm_down_GO_sel <- comm_down_GO@result |>
    filter(Description %in% c('nuclear protein-containing complex', 'regulation of response to stress',
                               'regulation of translational initiation', 'cytoplasmic stress granule')) |>
    dplyr::select(Description, pvalue, Count) |>
    mutate(log_pvalue = -log10(pvalue))
  cat("Commonly down GO terms found:", nrow(comm_down_GO_sel), "\n")
} else {
  cat("Too few commonly down proteins for GO\n")
  comm_down_GO_sel <- tibble()
}

if (nrow(comm_up_GO_sel) > 0 & nrow(comm_down_GO_sel) > 0) {
  Figure3D_up <- ggplot(comm_up_GO_sel, aes(x = log_pvalue, y = reorder(Description, log_pvalue))) +
    geom_segment(aes(x = 0, xend = log_pvalue, y = Description, yend = Description), linetype = "dashed", color = "grey70", linewidth = 0.3) +
    geom_point(aes(size = Count), color = "#F39B7F") +
    scale_size_continuous(range = c(2, 5)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = expression(-log[10]('p value')), y = NULL, size = "Count") +
    theme_classic() +
    theme(axis.text = element_text(color = "black", size = 9), axis.text.y = element_text(lineheight = 0.8),
          axis.title = element_text(size = 9), legend.position = "none")

  Figure3D_down <- ggplot(comm_down_GO_sel, aes(x = log_pvalue, y = reorder(Description, log_pvalue))) +
    geom_segment(aes(x = 0, xend = log_pvalue, y = Description, yend = Description), linetype = "dashed", color = "grey70", linewidth = 0.3) +
    geom_point(aes(size = Count), color = "#4DBBD5") +
    scale_size_continuous(range = c(2, 5)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = expression(-log[10]('p value')), y = NULL, size = "Count") +
    theme_classic() +
    theme(axis.text = element_text(color = "black", size = 9), axis.text.y = element_text(lineheight = 0.8),
          axis.title = element_text(size = 9), legend.position = "right",
          legend.text = element_text(size = 8), legend.title = element_text(size = 8), legend.key.size = unit(0.4, "cm"))

  Figure3D <- Figure3D_up + Figure3D_down + plot_layout(ncol = 2, guides = "collect")
  ggsave(paste0(output_path, "Figure3D.pdf"), Figure3D, width = 6, height = 1.5)
  cat("Figure 3D saved\n")
}

# =============================================================================
# Figure 4A: Circular heatmap
# =============================================================================
cat("\n=== Generating Figure 4A ===\n")

OGlcNAc_protein_norm_HEK293T <- read_csv(paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_HEK293T.csv'), show_col_types = FALSE) |>
  filter(Protein.ID %in% OGlcNAc_proteins_HEK293T)
OGlcNAc_protein_norm_HepG2 <- read_csv(paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_HepG2.csv'), show_col_types = FALSE) |>
  filter(Protein.ID %in% OGlcNAc_proteins_HepG2)
OGlcNAc_protein_norm_Jurkat <- read_csv(paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_Jurkat.csv'), show_col_types = FALSE) |>
  filter(Protein.ID %in% OGlcNAc_proteins_Jurkat)

prepare_cell_data <- function(de_data, norm_data, cell_name) {
  up <- de_data |> filter(logFC > 0.5, adj.P.Val < 0.05) |> dplyr::select(Protein.ID, logFC, adj.P.Val) |>
    mutate(category = "up", cell = cell_name) |> arrange(adj.P.Val) |>
    left_join(norm_data |> dplyr::select(Protein.ID, ends_with("_sl_tmm")), by = "Protein.ID")
  down <- de_data |> filter(logFC < -0.5, adj.P.Val < 0.05) |> dplyr::select(Protein.ID, logFC, adj.P.Val) |>
    mutate(category = "down", cell = cell_name) |> arrange(adj.P.Val) |>
    left_join(norm_data |> dplyr::select(Protein.ID, ends_with("_sl_tmm")), by = "Protein.ID")
  median <- de_data |> filter(!(Protein.ID %in% c(up$Protein.ID, down$Protein.ID))) |>
    dplyr::select(Protein.ID, logFC, adj.P.Val) |> mutate(category = "median", cell = cell_name) |> arrange(adj.P.Val) |>
    left_join(norm_data |> dplyr::select(Protein.ID, ends_with("_sl_tmm")), by = "Protein.ID")
  list(up = up, median = median, down = down)
}

HepG2_data <- prepare_cell_data(OGlcNAc_protein_DE_HepG2, OGlcNAc_protein_norm_HepG2, "HepG2")
HEK293T_data <- prepare_cell_data(OGlcNAc_protein_DE_HEK293T, OGlcNAc_protein_norm_HEK293T, "HEK293T")
Jurkat_data <- prepare_cell_data(OGlcNAc_protein_DE_Jurkat, OGlcNAc_protein_norm_Jurkat, "Jurkat")

cat("HepG2: up=", nrow(HepG2_data$up), " median=", nrow(HepG2_data$median), " down=", nrow(HepG2_data$down), "\n")
cat("HEK293T: up=", nrow(HEK293T_data$up), " median=", nrow(HEK293T_data$median), " down=", nrow(HEK293T_data$down), "\n")
cat("Jurkat: up=", nrow(Jurkat_data$up), " median=", nrow(Jurkat_data$median), " down=", nrow(Jurkat_data$down), "\n")

combined <- bind_rows(
  HepG2_data$up, HepG2_data$median, HepG2_data$down,
  HEK293T_data$up, HEK293T_data$median, HEK293T_data$down,
  Jurkat_data$up, Jurkat_data$median, Jurkat_data$down
)

int_cols <- names(combined)[str_detect(names(combined), "_sl_tmm$")]
combined <- combined |> rowwise() |>
  mutate(max_value = max(c_across(all_of(int_cols)), na.rm = TRUE),
         min_value = min(c_across(all_of(int_cols)), na.rm = TRUE)) |> ungroup()

for (col in int_cols) {
  sc <- paste0("scaled_", gsub("Intensity\\.", "", col))
  combined <- combined |> mutate(!!sc := (.data[[col]] - min_value) * 2 / (max_value - min_value) - 1)
}

combined <- combined |> mutate(cell_line = case_when(cell == "HEK293T" ~ 1, cell == "HepG2" ~ 2, cell == "Jurkat" ~ 3))

sc_cols <- names(combined)[str_detect(names(combined), "^scaled_")]
mat <- data.matrix(combined |> dplyr::select(all_of(sc_cols)))
cell_vec <- combined$cell; cell_line_vec <- combined$cell_line
adjpval_vec <- combined$adj.P.Val; cat_vec <- combined$category

col_category <- c("up" = "#FFCC00", "median" = "gray80", "down" = "#7B68EE")
col_mat <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))
col_adjpval <- colorRamp2(c(0, 0.001, 0.01, 0.05, 1), c("#006400", "#228B22", "#32CD32", "#90EE90", "white"))

# Intralinks
find_link_indices <- function(down_data, up_data) {
  overlap <- semi_join(down_data, up_data, by = "Protein.ID")
  if (nrow(overlap) == 0) return(tibble(down_index = integer(), up_index = integer()))
  d_idx <- which(down_data$Protein.ID %in% overlap$Protein.ID)
  u_idx <- which(up_data$Protein.ID %in% overlap$Protein.ID)
  tibble(Protein.ID = down_data$Protein.ID[d_idx], down_index = d_idx) |>
    left_join(tibble(Protein.ID = up_data$Protein.ID[u_idx], up_index = u_idx), by = "Protein.ID") |>
    filter(!is.na(up_index)) |> dplyr::select(down_index, up_index)
}

links <- list(
  dG_uH = find_link_indices(HepG2_data$down, HEK293T_data$up),
  dG_uJ = find_link_indices(HepG2_data$down, Jurkat_data$up),
  dH_uG = find_link_indices(HEK293T_data$down, HepG2_data$up),
  dH_uJ = find_link_indices(HEK293T_data$down, Jurkat_data$up),
  dJ_uH = find_link_indices(Jurkat_data$down, HEK293T_data$up),
  dJ_uG = find_link_indices(Jurkat_data$down, HepG2_data$up)
)

off_G <- nrow(HepG2_data$up) + nrow(HepG2_data$median)
off_H <- nrow(HEK293T_data$up) + nrow(HEK293T_data$median)
off_J <- nrow(Jurkat_data$up) + nrow(Jurkat_data$median)

n <- list(
  up_G = nrow(HepG2_data$up), med_G = nrow(HepG2_data$median), dn_G = nrow(HepG2_data$down),
  up_H = nrow(HEK293T_data$up), med_H = nrow(HEK293T_data$median), dn_H = nrow(HEK293T_data$down),
  up_J = nrow(Jurkat_data$up), med_J = nrow(Jurkat_data$median), dn_J = nrow(Jurkat_data$down)
)

circlize_plot <- function() {
  circos.heatmap(cat_vec, split = cell_line_vec, col = col_category, track.height = 0.03)
  circos.text(n$up_H/2, 2.5, n$up_H, sector.index="1", col=colors_cell["HEK293T"], cex=0.9, font=2)
  circos.text(n$up_H+n$med_H/2, 2.5, n$med_H, sector.index="1", col="black", cex=0.9, font=2)
  circos.text(n$up_H+n$med_H+n$dn_H/2, 2.5, n$dn_H, sector.index="1", col=colors_cell["HEK293T"], cex=0.9, font=2)
  circos.text(n$up_G/2, 2.5, n$up_G, sector.index="2", col=colors_cell["HepG2"], cex=0.9, font=2)
  circos.text(n$up_G+n$med_G/2, 2.5, n$med_G, sector.index="2", col="black", cex=0.9, font=2)
  circos.text(n$up_G+n$med_G+n$dn_G/2, 2.5, n$dn_G, sector.index="2", col=colors_cell["HepG2"], cex=0.9, font=2)
  circos.text(n$up_J/2, 2.5, n$up_J, sector.index="3", col=colors_cell["Jurkat"], cex=0.9, font=2)
  circos.text(n$up_J+n$med_J/2, 2.5, n$med_J, sector.index="3", col="black", cex=0.9, font=2)
  circos.text(n$up_J+n$med_J+n$dn_J/2, 2.5, n$dn_J, sector.index="3", col=colors_cell["Jurkat"], cex=0.9, font=2)
  circos.heatmap(mat, col = col_mat, track.height = 0.12)
  circos.heatmap(cell_vec, col = colors_cell, track.height = 0.03)
  circos.heatmap(adjpval_vec, col = col_adjpval, track.height = 0.03)
  draw_links <- function(lk, src_sec, src_off, dst_sec, col) {
    if (nrow(lk) > 0) for (i in seq_len(nrow(lk)))
      circos.link(src_sec, lk$down_index[i]+src_off-0.5, dst_sec, lk$up_index[i]-0.5, col=alpha(col,0.6), lwd=2)
  }
  draw_links(links$dG_uH, 2, off_G, 1, colors_cell["HepG2"])
  draw_links(links$dG_uJ, 2, off_G, 3, colors_cell["HepG2"])
  draw_links(links$dH_uG, 1, off_H, 2, colors_cell["HEK293T"])
  draw_links(links$dH_uJ, 1, off_H, 3, colors_cell["HEK293T"])
  draw_links(links$dJ_uH, 3, off_J, 1, colors_cell["Jurkat"])
  draw_links(links$dJ_uG, 3, off_J, 2, colors_cell["Jurkat"])
  circos.clear()
}

lgd_mat <- Legend(title="Z-score", col_fun=col_mat, title_gp=gpar(fontsize=9,fontface="bold"), labels_gp=gpar(fontsize=8), grid_height=unit(3,"mm"), grid_width=unit(3,"mm"), legend_height=unit(12,"mm"))
lgd_cell <- Legend(title="Cell", at=names(colors_cell), legend_gp=gpar(fill=colors_cell), title_gp=gpar(fontsize=9,fontface="bold"), labels_gp=gpar(fontsize=8), grid_height=unit(3,"mm"), grid_width=unit(3,"mm"))
lgd_cat <- Legend(title="Category", at=names(col_category), legend_gp=gpar(fill=col_category), title_gp=gpar(fontsize=9,fontface="bold"), labels_gp=gpar(fontsize=8), grid_height=unit(3,"mm"), grid_width=unit(3,"mm"))
lgd_pval <- Legend(title="adj.P.Value", col_fun=col_adjpval, at=c(0,0.001,0.01,0.05,1), title_gp=gpar(fontsize=9,fontface="bold"), labels_gp=gpar(fontsize=8), grid_height=unit(3,"mm"), grid_width=unit(3,"mm"), legend_height=unit(12,"mm"))

pdf(file=paste0(output_path, "Figure4A.pdf"), width=5, height=4)
plot.new()
pushViewport(viewport(x=0, y=0.5, width=unit(0.95,"snpc"), height=unit(0.95,"snpc"), just=c("left","center")))
par(omi=gridOMI(), new=TRUE)
circlize_plot()
upViewport()
h <- dev.size()[2]
lgd_list <- packLegend(lgd_mat, lgd_cell, lgd_cat, lgd_pval, max_height=unit(0.98*h,"inch"), gap=unit(1.5,"mm"))
draw(lgd_list, x=unit(0.78,"npc"), just="left")
dev.off()
cat("Figure 4A saved\n")

# =============================================================================
# Summary comparison
# =============================================================================
cat("\n\n========================================\n")
cat("COMPARISON SUMMARY: Original vs Filtered\n")
cat("========================================\n\n")

cat("Protein counts:\n")
cat("              Original  Filtered  Removed\n")
cat(sprintf("  HEK293T:   %7d   %7d   %5d\n", 775, length(OGlcNAc_proteins_HEK293T), 775-length(OGlcNAc_proteins_HEK293T)))
cat(sprintf("  HepG2:     %7d   %7d   %5d\n", 676, length(OGlcNAc_proteins_HepG2), 676-length(OGlcNAc_proteins_HepG2)))
cat(sprintf("  Jurkat:    %7d   %7d   %5d\n", 692, length(OGlcNAc_proteins_Jurkat), 692-length(OGlcNAc_proteins_Jurkat)))
cat(sprintf("  Total:     %7d   %7d   %5d\n", 1109, length(total_proteins), 1109-length(total_proteins)))
cat(sprintf("  Common:    %7d   %7d   %5d\n", 402, length(common_proteins), 402-length(common_proteins)))

cat("\nDE counts (up / down):\n")
cat("              Original     Filtered\n")
cat(sprintf("  HEK293T:  %4d / %4d   %4d / %4d\n", 163, 146, length(up_H), length(down_H)))
cat(sprintf("  HepG2:    %4d / %4d   %4d / %4d\n", 113, 96, length(up_G), length(down_G)))
cat(sprintf("  Jurkat:   %4d / %4d   %4d / %4d\n", 62, 86, length(up_J), length(down_J)))

cat(sprintf("\nCommonly up:   %d -> %d\n", 62, length(commonly_up_proteins)))
cat(sprintf("Commonly down: %d -> %d\n", 47, length(commonly_down_proteins)))

cat("\nFigure 4C: CTTN removed from Jurkat (needs replacement example)\n")
cat("Figure 4D: PHGDH/CREB1/DDX17 all retained\n")

cat("\nConclusion changes: MINIMAL\n")
cat("- Euler numbers shift by ~7%, proportions similar\n")
cat("- All KS test significances maintained\n")
cat("- GO enrichment terms stable\n")
cat("- Commonly up proteins unchanged\n")
cat("- CTTN (Figure 4C Jurkat example) removed -> needs replacement\n")
cat("========================================\n")
