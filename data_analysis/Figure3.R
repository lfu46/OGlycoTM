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
    axis.title = element_text(size = 7),
    axis.text.x = element_text(size = 7, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 7, color = "black"),
    legend.position = "none"
  )

print(Figure3A)

ggsave(
  filename = paste0(figure_file_path, "Figure3/Figure3A.pdf"),
  plot = Figure3A,
  width = 1.5, height = 1.5, units = "in"
)

cat("\nFigure 3A saved to:", figure_file_path, "Figure3/\n")

# Figure 3B ---------------------------------------------------------------
# (Add later)

# Figure 3C ---------------------------------------------------------------
# Violin boxplot comparing logFC of O-GlcNAc proteins in specific functional categories
# Categories: DNA binding, RNA binding, Transcription
# Statistical test: Kolmogorov-Smirnov test (compares distribution shapes)

# Read GO enrichment results for common O-GlcNAc proteins
common_OGlcNAc_GO <- read_csv(
  paste0(source_file_path, 'enrichment/Figure2C_common_OGlcNAc_GO.csv')
)

# Define GO terms for each functional category
DNA_binding_terms <- c(
  "chromatin DNA binding",
  "damaged DNA binding",
  "telomeric DNA binding",
  "core promoter sequence-specific DNA binding",
  "single-stranded DNA binding"
)

RNA_binding_terms <- c(
  "mRNA binding",
  "single-stranded RNA binding",
  "double-stranded RNA binding",
  "pre-mRNA binding",
  "telomerase RNA binding",
  "snRNA binding"
)

# Extract all transcription-related terms using pattern matching
transcription_terms <- common_OGlcNAc_GO %>%
  filter(str_detect(Description, regex("transcription", ignore_case = TRUE))) %>%
  pull(Description)

# Extract proteins for each category
extract_proteins <- function(go_df, terms) {
  go_df %>%
    filter(Description %in% terms) %>%
    pull(geneID) %>%
    paste(collapse = "/") %>%
    str_split("/") %>%
    unlist() %>%
    unique()
}

DNA_binding_proteins <- extract_proteins(common_OGlcNAc_GO, DNA_binding_terms)
RNA_binding_proteins <- extract_proteins(common_OGlcNAc_GO, RNA_binding_terms)
transcription_proteins <- extract_proteins(common_OGlcNAc_GO, transcription_terms)

cat("Figure 3C - DNA binding proteins:", length(DNA_binding_proteins), "\n")
cat("Figure 3C - RNA binding proteins:", length(RNA_binding_proteins), "\n")
cat("Figure 3C - Transcription proteins:", length(transcription_proteins), "\n")

# Create combined logFC data for each category
create_category_df <- function(proteins, category_name) {
  bind_rows(
    OGlcNAc_protein_DE_HEK293T %>%
      filter(Protein.ID %in% proteins) %>%
      select(Protein.ID, logFC) %>%
      mutate(CellType = "HEK293T"),
    OGlcNAc_protein_DE_HepG2 %>%
      filter(Protein.ID %in% proteins) %>%
      select(Protein.ID, logFC) %>%
      mutate(CellType = "HepG2"),
    OGlcNAc_protein_DE_Jurkat %>%
      filter(Protein.ID %in% proteins) %>%
      select(Protein.ID, logFC) %>%
      mutate(CellType = "Jurkat")
  ) %>%
    mutate(Category = category_name)
}

Figure3C_df <- bind_rows(
  create_category_df(DNA_binding_proteins, "DNA binding"),
  create_category_df(RNA_binding_proteins, "RNA binding"),
  create_category_df(transcription_proteins, "Transcription")
)

# Set factor levels
Figure3C_df <- Figure3C_df %>%
  mutate(
    CellType = factor(CellType, levels = c("HEK293T", "HepG2", "Jurkat")),
    Category = factor(Category, levels = c("DNA binding", "RNA binding", "Transcription"))
  )

# Perform KS tests for each category (between cell types)
perform_ks_tests <- function(df, category_name) {
  df_cat <- df %>% filter(Category == category_name)

  HEK293T_logFC <- df_cat %>% filter(CellType == "HEK293T") %>% pull(logFC)
  HepG2_logFC <- df_cat %>% filter(CellType == "HepG2") %>% pull(logFC)
  Jurkat_logFC <- df_cat %>% filter(CellType == "Jurkat") %>% pull(logFC)

  tribble(
    ~Category, ~group1, ~group2, ~p,
    category_name, "HEK293T", "HepG2",
    if(length(HEK293T_logFC) > 0 & length(HepG2_logFC) > 0) ks.test(HEK293T_logFC, HepG2_logFC)$p.value else NA,
    category_name, "HEK293T", "Jurkat",
    if(length(HEK293T_logFC) > 0 & length(Jurkat_logFC) > 0) ks.test(HEK293T_logFC, Jurkat_logFC)$p.value else NA,
    category_name, "HepG2", "Jurkat",
    if(length(HepG2_logFC) > 0 & length(Jurkat_logFC) > 0) ks.test(HepG2_logFC, Jurkat_logFC)$p.value else NA
  )
}

Figure3C_ks_results <- bind_rows(
  perform_ks_tests(Figure3C_df, "DNA binding"),
  perform_ks_tests(Figure3C_df, "RNA binding"),
  perform_ks_tests(Figure3C_df, "Transcription")
) %>%
  filter(!is.na(p)) %>%
  add_significance("p") %>%
  filter(p.signif != "ns")  # Remove non-significant comparisons

cat("\nFigure 3C - KS Test Results:\n")
print(Figure3C_ks_results)

# Add y positions for significance bars
Figure3C_ks_results <- Figure3C_ks_results %>%
  mutate(
    y.position = case_when(
      group1 == "HEK293T" & group2 == "HepG2" ~ 1.0,
      group1 == "HEK293T" & group2 == "Jurkat" ~ 1.4,
      group1 == "HepG2" & group2 == "Jurkat" ~ 1.2,
      TRUE ~ 1.0
    )
  )

# Create faceted violin boxplot
Figure3C <- Figure3C_df %>%
  ggplot(aes(x = CellType, y = logFC)) +
  geom_violin(aes(fill = CellType), color = "transparent") +
  geom_boxplot(color = "black", outliers = FALSE, width = 0.2, linewidth = 0.3) +
  facet_wrap(~ Category, nrow = 1) +
  scale_fill_manual(values = colors_cell) +
  labs(
    x = "",
    y = expression(log[2]*"(Tuni/Ctrl)")
  ) +
  stat_pvalue_manual(
    data = Figure3C_ks_results,
    label = "p.signif",
    tip.length = 0,
    size = 3
  ) +
  coord_cartesian(ylim = c(-1.5, 1.7)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 7),
    axis.text.x = element_text(size = 6, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 7, color = "black"),
    strip.text = element_text(size = 6),
    strip.background = element_blank(),
    legend.position = "none"
  )

print(Figure3C)

ggsave(
  filename = paste0(figure_file_path, "Figure3/Figure3C.pdf"),
  plot = Figure3C,
  width = 2, height = 1.5, units = "in"
)

cat("\nFigure 3C saved to:", figure_file_path, "Figure3/\n")

# Figure 3D ---------------------------------------------------------------
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

cat("Figure 3D - Number of overlapping proteins:", length(overlapping_proteins), "\n")

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
Figure3D <- ggplot() +
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
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7, color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.15, "cm")
  )

print(Figure3D)

ggsave(
  filename = paste0(figure_file_path, "Figure3/Figure3D.pdf"),
  plot = Figure3D,
  width = 2.5, height = 1.5, units = "in"
)

cat("\nFigure 3D saved to:", figure_file_path, "Figure3/\n")

# Figure 3E ---------------------------------------------------------------
# GO enrichment analysis for commonly up/downregulated O-GlcNAc proteins
# Comparing two backgrounds: overlapping proteins (402) vs total O-GlcNAc (1109)

library(clusterProfiler)
library(org.Hs.eg.db)

# Get overlapping proteins (if not already created in Figure 3D)
if (!exists("overlapping_proteins")) {
  overlapping_proteins <- Reduce(intersect, list(
    OGlcNAc_protein_DE_HEK293T$Protein.ID,
    OGlcNAc_protein_DE_HepG2$Protein.ID,
    OGlcNAc_protein_DE_Jurkat$Protein.ID
  ))
}

# Get protein lists for enrichment
commonly_up_proteins <- OGlcNAc_protein_commonly_up$Protein.ID
commonly_down_proteins <- OGlcNAc_protein_commonly_down$Protein.ID

cat("Figure 3E - Commonly upregulated proteins:", length(commonly_up_proteins), "\n")
cat("Figure 3E - Commonly downregulated proteins:", length(commonly_down_proteins), "\n")

# Background 1: Overlapping O-GlcNAc proteins (402)
background_overlapping <- overlapping_proteins
cat("Background (overlapping):", length(background_overlapping), "proteins\n")

# Background 2: Total quantified O-GlcNAc proteins (union of all 3 cell types)
background_total <- Reduce(union, list(
  OGlcNAc_protein_DE_HEK293T$Protein.ID,
  OGlcNAc_protein_DE_HepG2$Protein.ID,
  OGlcNAc_protein_DE_Jurkat$Protein.ID
))
cat("Background (total quantified):", length(background_total), "proteins\n")

# --- Enrichment with overlapping background ---

# GO enrichment for commonly upregulated proteins
commonly_up_GO_overlapping <- enrichGO(
  gene = commonly_up_proteins,
  OrgDb = org.Hs.eg.db,
  universe = background_overlapping,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  commonly_up_GO_overlapping@result,
  paste0(source_file_path, 'enrichment/Figure3E_commonly_up_GO_overlapping.csv')
)

# GO enrichment for commonly downregulated proteins
commonly_down_GO_overlapping <- enrichGO(
  gene = commonly_down_proteins,
  OrgDb = org.Hs.eg.db,
  universe = background_overlapping,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  commonly_down_GO_overlapping@result,
  paste0(source_file_path, 'enrichment/Figure3E_commonly_down_GO_overlapping.csv')
)

# --- Enrichment with total O-GlcNAc background ---

# GO enrichment for commonly upregulated proteins
commonly_up_GO_total <- enrichGO(
  gene = commonly_up_proteins,
  OrgDb = org.Hs.eg.db,
  universe = background_total,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  commonly_up_GO_total@result,
  paste0(source_file_path, 'enrichment/Figure3E_commonly_up_GO_total.csv')
)

# GO enrichment for commonly downregulated proteins
commonly_down_GO_total <- enrichGO(
  gene = commonly_down_proteins,
  OrgDb = org.Hs.eg.db,
  universe = background_total,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  commonly_down_GO_total@result,
  paste0(source_file_path, 'enrichment/Figure3E_commonly_down_GO_total.csv')
)

cat("\nEnrichment results saved to:", source_file_path, "enrichment/\n")

# --- Create dot plot (using total O-GlcNAc background results) ---

# Select GO terms for commonly upregulated proteins (manually select 3-4 terms)
commonly_up_GO_selected <- read_csv(
  paste0(source_file_path, 'enrichment/Figure3E_commonly_up_GO_total.csv')
) |>
  filter(
    Description %in% c(
      # TODO: Add 3-4 GO terms for commonly upregulated proteins
      'organophosphate metabolic process',
      'regulation of translation',
      'response to glucose',
      'positive regulation of RNA splicing'
    )
  ) |>
  dplyr::select(Description, ONTOLOGY, pvalue, Count) |>
  mutate(
    log_pvalue = -log10(pvalue),
    category = "Commonly Up"
  )

# Select GO terms for commonly downregulated proteins (manually select 3-4 terms)
commonly_down_GO_selected <- read_csv(
  paste0(source_file_path, 'enrichment/Figure3E_commonly_down_GO_total.csv')
) |>
  filter(
    Description %in% c(
      # TODO: Add 3-4 GO terms for commonly downregulated proteins
      'nuclear protein-containing complex',
      'regulation of response to stress',
      'regulation of translational initiation',
      'cytoplasmic stress granule'
    )
  ) |>
  dplyr::select(Description, ONTOLOGY, pvalue, Count) |>
  mutate(
    log_pvalue = -log10(pvalue),
    category = "Commonly Down"
  )

# Combine for plotting
Figure3E_GO_df <- bind_rows(
  commonly_up_GO_selected,
  commonly_down_GO_selected
) |>
  mutate(
    category = factor(category, levels = c("Commonly Up", "Commonly Down"))
  )

# Define colors
Figure3E_colors <- c(
  "Commonly Up" = "#F39B7F",
  "Commonly Down" = "#4DBBD5"
)

# Create dot plot
Figure3E <- ggplot(Figure3E_GO_df, aes(x = log_pvalue, y = Description)) +
  geom_segment(
    aes(x = 0, xend = Inf, y = Description, yend = Description),
    linetype = "dashed",
    color = "grey70",
    linewidth = 0.3
  ) +
  geom_point(aes(size = Count, color = category)) +
  facet_grid(rows = vars(category), scales = "free_y", space = "free_y") +
  scale_color_manual(values = Figure3E_colors) +
  scale_size_continuous(range = c(1, 4)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = expression(-log[10]('p value')),
    y = NULL,
    size = "Count"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 7),
    axis.text.y = element_text(lineheight = 0.6),
    axis.title = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(size = 5),
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    legend.key.size = unit(0.3, "cm")
  ) +
  guides(color = "none")

print(Figure3E)

ggsave(
  filename = paste0(figure_file_path, "Figure3/Figure3E.pdf"),
  plot = Figure3E,
  width = 3, height = 1.5, units = "in"
)

cat("\nFigure 3E saved to:", figure_file_path, "Figure3/\n")

