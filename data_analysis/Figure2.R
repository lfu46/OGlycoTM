# import packages
library(tidyverse)

# Figure 2. Identification of O-GlcNAc and O-GalNAc sites and proteins in HEK293T, HepG2 and Jurkat cells


# Figure 2A ---------------------------------------------------------------
# Identification of O-GlcNAc proteins in three types of cells
# Stacked bar: common (grey) at bottom, cell-type specific (colored) on top

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

# Calculate common O-GlcNAc proteins (intersection of all 3 cell types)
common_OGlcNAc_proteins <- intersect(
  intersect(OGlcNAc_proteins_HEK293T, OGlcNAc_proteins_HepG2),
  OGlcNAc_proteins_Jurkat
)
common_count <- length(common_OGlcNAc_proteins)

# Calculate total and non-common counts for each cell type
HEK293T_total <- length(OGlcNAc_proteins_HEK293T)
HepG2_total <- length(OGlcNAc_proteins_HepG2)
Jurkat_total <- length(OGlcNAc_proteins_Jurkat)

HEK293T_noncommon <- HEK293T_total - common_count
HepG2_noncommon <- HepG2_total - common_count
Jurkat_noncommon <- Jurkat_total - common_count

# Create data frame for stacked bar plot
# Common (grey) at bottom, cell-type specific (colored) on top
Figure2A_df <- data.frame(
  cell_type = factor(
    rep(c("HEK293T", "HepG2", "Jurkat"), each = 2),
    levels = c("HEK293T", "HepG2", "Jurkat")
  ),
  category = factor(
    rep(c("Common", "Cell-specific"), 3),
    levels = c("Common", "Cell-specific")
  ),
  protein_count = c(
    common_count, HEK293T_noncommon,
    common_count, HepG2_noncommon,
    common_count, Jurkat_noncommon
  ),
  fill_color = factor(
    c(
      "Common", "HEK293T",
      "Common", "HepG2",
      "Common", "Jurkat"
    ),
    levels = c("Common", "HEK293T", "HepG2", "Jurkat")
  )
)

# Define colors: grey for common, cell-type colors for specific
Figure2A_colors <- c(
  "Common" = "grey70",
  "HEK293T" = "#4DBBD5",
  "HepG2" = "#F39B7F",
  "Jurkat" = "#00A087"
)

# Calculate total counts for labels
Figure2A_totals <- data.frame(
  cell_type = factor(c("HEK293T", "HepG2", "Jurkat"), levels = c("HEK293T", "HepG2", "Jurkat")),
  total = c(HEK293T_total, HepG2_total, Jurkat_total)
)

# Stacked barplot
Figure2A <- ggplot(Figure2A_df, aes(x = cell_type, y = protein_count, fill = fill_color)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  geom_text(
    data = Figure2A_totals,
    aes(x = cell_type, y = total, label = total, fill = NULL),
    vjust = -0.2,
    size = 2,
    color = "black"
  ) +
  scale_fill_manual(
    values = Figure2A_colors,
    breaks = c("Common", "HEK293T", "HepG2", "Jurkat"),
    labels = c("Common", "HEK293T", "HepG2", "Jurkat")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = NULL,
    y = "No. of O-GlcNAc proteins",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    title = element_text(family = "Helvetica", color = "black", size = 6),
    text = element_text(family = "Helvetica", color = "black"),
    axis.text = element_text(color = "black", size = 6),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.key.size = unit(0.2, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6)
  )

ggsave(
  filename = paste0(figure_file_path, "Figure2/Figure2A.pdf"),
  plot = Figure2A,
  width = 2, height = 1.5, units = "in"
)


# Figure 2B ---------------------------------------------------------------
# Localized O-GlcNAc sites in three types of cells

# generate site index for localized O-GlcNAc and O-GalNAc
# HEK293T
OGlyco_site_HEK293T <- OGlyco_HEK293T_bonafide |>
  filter(
    Confidence.Level %in% c('Level1', 'Level1b')
  ) |>
  mutate(
    peptide_site = as.numeric(str_extract(Site.Probabilities, "(?<=\\[)\\d+")),
    modified_residue = substr(Peptide, peptide_site, peptide_site),
    site_number = Protein.Start + peptide_site - 1,
    site_index = paste0(Protein.ID, "_", modified_residue, site_number),
    glycan_type = case_when(
      Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859') ~ "O-GlcNAc",
      Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1) % 326.1339', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968') ~ "O-GalNAc"
    )
  )

# HepG2
OGlyco_site_HepG2 <- OGlyco_HepG2_bonafide |>
  filter(
    Confidence.Level %in% c('Level1', 'Level1b')
  ) |>
  mutate(
    peptide_site = as.numeric(str_extract(Site.Probabilities, "(?<=\\[)\\d+")),
    modified_residue = substr(Peptide, peptide_site, peptide_site),
    site_number = Protein.Start + peptide_site - 1,
    site_index = paste0(Protein.ID, "_", modified_residue, site_number),
    glycan_type = case_when(
      Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859') ~ "O-GlcNAc",
      Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1) % 326.1339', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968') ~ "O-GalNAc"
    )
  )

# Jurkat
OGlyco_site_Jurkat <- OGlyco_Jurkat_bonafide |>
  filter(
    Confidence.Level %in% c('Level1', 'Level1b')
  ) |>
  mutate(
    peptide_site = as.numeric(str_extract(Site.Probabilities, "(?<=\\[)\\d+")),
    modified_residue = substr(Peptide, peptide_site, peptide_site),
    site_number = Protein.Start + peptide_site - 1,
    site_index = paste0(Protein.ID, "_", modified_residue, site_number),
    glycan_type = case_when(
      Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859') ~ "O-GlcNAc",
      Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1) % 326.1339', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968') ~ "O-GalNAc"
    )
  )

# save to csv
write_csv(OGlyco_site_HEK293T, paste0(source_file_path, "site/OGlyco_site_HEK293T.csv"))
write_csv(OGlyco_site_HepG2, paste0(source_file_path, "site/OGlyco_site_HepG2.csv"))
write_csv(OGlyco_site_Jurkat, paste0(source_file_path, "site/OGlyco_site_Jurkat.csv"))

# Venn diagram for O-GlcNAc sites across three cell types
library(eulerr)

# pull unique O-GlcNAc site index from each cell type
OGlcNAc_site_HEK293T <- OGlyco_site_HEK293T |>
  filter(glycan_type == "O-GlcNAc") |>
  pull(site_index) |>
  unique()

OGlcNAc_site_HepG2 <- OGlyco_site_HepG2 |>
  filter(glycan_type == "O-GlcNAc") |>
  pull(site_index) |>
  unique()

OGlcNAc_site_Jurkat <- OGlyco_site_Jurkat |>
  filter(glycan_type == "O-GlcNAc") |>
  pull(site_index) |>
  unique()

# create list for Venn diagram with counts in labels
OGlcNAc_site_list <- list(
  OGlcNAc_site_HEK293T,
  OGlcNAc_site_HepG2,
  OGlcNAc_site_Jurkat
)
names(OGlcNAc_site_list) <- c(
  paste0("HEK293T (", length(OGlcNAc_site_HEK293T), ")"),
  paste0("HepG2 (", length(OGlcNAc_site_HepG2), ")"),
  paste0("Jurkat (", length(OGlcNAc_site_Jurkat), ")")
)

# create euler object (proportional)
OGlcNAc_euler <- euler(OGlcNAc_site_list)

# Venn diagram
# Create color vector with updated names
Figure2B_colors <- c(
  colors_cell["HEK293T"],
  colors_cell["HepG2"],
  colors_cell["Jurkat"]
)
names(Figure2B_colors) <- names(OGlcNAc_site_list)

Figure2B <- plot(
  OGlcNAc_euler,
  fills = list(fill = Figure2B_colors, alpha = 0.5),
  edges = list(col = "white", lwd = 2),
  labels = list(font = 1, cex = 0.5),
  quantities = list(font = 1, cex = 0.5)
)

# save plot
pdf(paste0(figure_file_path, "Figure2/Figure2B.pdf"), width = 2, height = 1.5)
print(Figure2B)
dev.off()


# Figure 2C ---------------------------------------------------------------
# Donut chart showing distribution of O-GlcNAc proteins across cell types
# Categories: Common (all 3), Unique (each cell type), Shared (exactly 2 cell types)
# Functional enrichment analysis for common and unique proteins

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

# Calculate total O-GlcNAc proteins (union of all 3 cell types)
total_OGlcNAc_proteins <- union(
  union(OGlcNAc_proteins_HEK293T, OGlcNAc_proteins_HepG2),
  OGlcNAc_proteins_Jurkat
)
total_OGlcNAc_count <- length(total_OGlcNAc_proteins)

# Calculate common O-GlcNAc proteins (intersection of all 3 cell types)
common_OGlcNAc_proteins <- intersect(
  intersect(OGlcNAc_proteins_HEK293T, OGlcNAc_proteins_HepG2),
  OGlcNAc_proteins_Jurkat
)
common_OGlcNAc_count <- length(common_OGlcNAc_proteins)

# Calculate unique O-GlcNAc proteins for each cell type (only in that cell type)
unique_OGlcNAc_HEK293T <- setdiff(
  setdiff(OGlcNAc_proteins_HEK293T, OGlcNAc_proteins_HepG2),
  OGlcNAc_proteins_Jurkat
)
unique_OGlcNAc_HEK293T_count <- length(unique_OGlcNAc_HEK293T)

unique_OGlcNAc_HepG2 <- setdiff(
  setdiff(OGlcNAc_proteins_HepG2, OGlcNAc_proteins_HEK293T),
  OGlcNAc_proteins_Jurkat
)
unique_OGlcNAc_HepG2_count <- length(unique_OGlcNAc_HepG2)

unique_OGlcNAc_Jurkat <- setdiff(
  setdiff(OGlcNAc_proteins_Jurkat, OGlcNAc_proteins_HEK293T),
  OGlcNAc_proteins_HepG2
)
unique_OGlcNAc_Jurkat_count <- length(unique_OGlcNAc_Jurkat)

# Calculate shared O-GlcNAc proteins (exactly 2 cell types, not all 3)
# HEK293T & HepG2 only (not in Jurkat)
shared_HEK293T_HepG2 <- setdiff(
  intersect(OGlcNAc_proteins_HEK293T, OGlcNAc_proteins_HepG2),
  OGlcNAc_proteins_Jurkat
)

# HEK293T & Jurkat only (not in HepG2)
shared_HEK293T_Jurkat <- setdiff(
  intersect(OGlcNAc_proteins_HEK293T, OGlcNAc_proteins_Jurkat),
  OGlcNAc_proteins_HepG2
)

# HepG2 & Jurkat only (not in HEK293T)
shared_HepG2_Jurkat <- setdiff(
  intersect(OGlcNAc_proteins_HepG2, OGlcNAc_proteins_Jurkat),
  OGlcNAc_proteins_HEK293T
)

# Combined shared proteins (all proteins in exactly 2 cell types)
shared_OGlcNAc_proteins <- union(
  union(shared_HEK293T_HepG2, shared_HEK293T_Jurkat),
  shared_HepG2_Jurkat
)
shared_OGlcNAc_count <- length(shared_OGlcNAc_proteins)

# Verify total count
# total should equal: common + unique_HEK293T + unique_HepG2 + unique_Jurkat + shared
verify_total <- common_OGlcNAc_count + unique_OGlcNAc_HEK293T_count +
  unique_OGlcNAc_HepG2_count + unique_OGlcNAc_Jurkat_count + shared_OGlcNAc_count
cat("Total O-GlcNAc proteins:", total_OGlcNAc_count, "\n")
cat("Sum of categories:", verify_total, "\n")

# Save protein lists to CSV
write_csv(
  data.frame(Protein.ID = total_OGlcNAc_proteins),
  paste0(source_file_path, 'protein_lists/OGlcNAc_protein_total.csv')
)
write_csv(
  data.frame(Protein.ID = common_OGlcNAc_proteins),
  paste0(source_file_path, 'protein_lists/common_OGlcNAc_proteins.csv')
)
write_csv(
  data.frame(Protein.ID = unique_OGlcNAc_HEK293T),
  paste0(source_file_path, 'protein_lists/unique_OGlcNAc_HEK293T.csv')
)
write_csv(
  data.frame(Protein.ID = unique_OGlcNAc_HepG2),
  paste0(source_file_path, 'protein_lists/unique_OGlcNAc_HepG2.csv')
)
write_csv(
  data.frame(Protein.ID = unique_OGlcNAc_Jurkat),
  paste0(source_file_path, 'protein_lists/unique_OGlcNAc_Jurkat.csv')
)
write_csv(
  data.frame(Protein.ID = shared_OGlcNAc_proteins),
  paste0(source_file_path, 'protein_lists/shared_OGlcNAc_proteins.csv')
)

# Define colors for donut chart
Figure2C_colors <- c(
  "Common" = "#808080",
  "HEK293T" = "#4DBBD5",
  "HepG2" = "#F39B7F",
  "Jurkat" = "#00A087",
  "Shared" = "#B8A9C9"
)

# Create data for donut chart
# Order: HEK293T -> HepG2 -> Jurkat -> Shared -> Common (clockwise from 12 o'clock)
Figure2C_df <- data.frame(
  category = factor(
    c("HEK293T", "HepG2", "Jurkat", "Shared", "Common"),
    levels = c("HEK293T", "HepG2", "Jurkat", "Shared", "Common")
  ),
  count = c(
    unique_OGlcNAc_HEK293T_count,
    unique_OGlcNAc_HepG2_count,
    unique_OGlcNAc_Jurkat_count,
    shared_OGlcNAc_count,
    common_OGlcNAc_count
  )
)

# Calculate percentages and positions for labels
Figure2C_df <- Figure2C_df |>
  mutate(
    percentage = count / sum(count) * 100,
    ymax = cumsum(percentage),
    ymin = lag(ymax, default = 0),
    label_pos = (ymin + ymax) / 2
  )

# Donut chart
Figure2C_donut <- ggplot(Figure2C_df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.5, fill = category)) +
  geom_rect(color = "white", linewidth = 0.3) +
  geom_text(
    aes(x = 3.25, y = label_pos, label = count),
    size = 2,
    color = "black"
  ) +
  # Add "O-GlcNAc" label in center
  annotate(
    "text",
    x = 0, y = 0,
    label = "O-GlcNAc",
    size = 2.5,
    fontface = "bold"
  ) +
  coord_polar(theta = "y") +
  xlim(c(0, 4)) +
  # Remove "Shared" from legend but keep in plot
  scale_fill_manual(
    values = Figure2C_colors,
    breaks = c("Common", "HEK293T", "HepG2", "Jurkat")
  ) +
  labs(fill = NULL) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.2, "cm"),
    legend.text = element_text(family = "Helvetica", size = 5)
  )

# Functional enrichment analysis using clusterProfiler
library(clusterProfiler)
library(org.Hs.eg.db)
library(Biostrings)

# Get reviewed human proteome as background from UniProt FASTA file
human_proteome_fasta <- readAAStringSet(
  paste0(source_file_path, 'reference/uniprotkb_reviewed_true_AND_model_organ_2026_01_09.fasta')
)

# Extract UniProt IDs from FASTA headers (format: >sp|P12345|GENE_HUMAN ...)
human_proteome <- names(human_proteome_fasta) |>
  str_extract("(?<=\\|)[A-Z0-9]+(?=\\|)")

## GO enrichment for common O-GlcNAc proteins (background: human proteome)
common_OGlcNAc_GO <- enrichGO(

  gene = common_OGlcNAc_proteins,
  OrgDb = org.Hs.eg.db,
  universe = human_proteome,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  common_OGlcNAc_GO@result,
  paste0(source_file_path, 'enrichment/Figure2C_common_OGlcNAc_GO.csv')
)

## GO enrichment for unique HEK293T O-GlcNAc proteins (background: total O-GlcNAc)
unique_HEK293T_OGlcNAc_GO <- enrichGO(
  gene = unique_OGlcNAc_HEK293T,
  OrgDb = org.Hs.eg.db,
  universe = total_OGlcNAc_proteins,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  unique_HEK293T_OGlcNAc_GO@result,
  paste0(source_file_path, 'enrichment/Figure2C_unique_HEK293T_OGlcNAc_GO.csv')
)

## GO enrichment for unique HepG2 O-GlcNAc proteins (background: total O-GlcNAc)
unique_HepG2_OGlcNAc_GO <- enrichGO(
  gene = unique_OGlcNAc_HepG2,
  OrgDb = org.Hs.eg.db,
  universe = total_OGlcNAc_proteins,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  unique_HepG2_OGlcNAc_GO@result,
  paste0(source_file_path, 'enrichment/Figure2C_unique_HepG2_OGlcNAc_GO.csv')
)

## GO enrichment for unique Jurkat O-GlcNAc proteins (background: total O-GlcNAc)
unique_Jurkat_OGlcNAc_GO <- enrichGO(
  gene = unique_OGlcNAc_Jurkat,
  OrgDb = org.Hs.eg.db,
  universe = total_OGlcNAc_proteins,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  unique_Jurkat_OGlcNAc_GO@result,
  paste0(source_file_path, 'enrichment/Figure2C_unique_Jurkat_OGlcNAc_GO.csv')
)

# Import and select GO terms for dot plot (fill in terms after reviewing results)
# Common O-GlcNAc proteins (separate plot on left side)
# Select 3 terms per ontology (BP, MF, CC) = 9 terms total
common_OGlcNAc_GO_selected <- read_csv(
  paste0(source_file_path, 'enrichment/Figure2C_common_OGlcNAc_GO.csv')
) |>
  filter(
    Description %in% c(
      # TODO: Add 3 BP terms here
      'regulation of mRNA metabolic process',
      'nucleocytoplasmic transport',
      'nuclear transport',
      # TODO: Add 3 MF terms here
      'transcription coregulator activity',
      'DNA-binding transcription factor binding',
      'mRNA binding',
      # TODO: Add 3 CC terms here
      'cytoplasmic ribonucleoprotein granule',
      'histone acetyltransferase complex',
      'cytoplasmic stress granule'
    )
  ) |>
  dplyr::select(ONTOLOGY, Description, pvalue, Count) |>
  mutate(
    log_pvalue = -log10(pvalue),
    ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "MF", "CC")),
    Description = str_wrap(Description, width = 30)
  )

# Unique HEK293T O-GlcNAc proteins (select 3 terms)
unique_HEK293T_GO_selected <- read_csv(
  paste0(source_file_path, 'enrichment/Figure2C_unique_HEK293T_OGlcNAc_GO.csv')
) |>
  filter(
    Description %in% c(
      # TODO: Add 3 GO terms here for HEK293T
      'chaperone-mediated protein folding',
      'cell cycle process',
      'epithelial cell proliferation'
    )
  ) |>
  dplyr::select(Description, pvalue, Count) |>
  mutate(
    log_pvalue = -log10(pvalue),
    category = "HEK293T"
  )

# Unique HepG2 O-GlcNAc proteins (select 3 terms)
unique_HepG2_GO_selected <- read_csv(
  paste0(source_file_path, 'enrichment/Figure2C_unique_HepG2_OGlcNAc_GO.csv')
) |>
  filter(
    Description %in% c(
      # TODO: Add 3 GO terms here for HepG2
      'extracellular matrix',
      'cell periphery',
      'endoplasmic reticulum lumen'
    )
  ) |>
  dplyr::select(Description, pvalue, Count) |>
  mutate(
    log_pvalue = -log10(pvalue),
    category = "HepG2"
  )

# Unique Jurkat O-GlcNAc proteins (select 3 terms)
unique_Jurkat_GO_selected <- read_csv(
  paste0(source_file_path, 'enrichment/Figure2C_unique_Jurkat_OGlcNAc_GO.csv')
) |>
  filter(
    Description %in% c(
      # TODO: Add 3 GO terms here for Jurkat
      'leukocyte activation involved in immune response',
      'leukocyte cell-cell adhesion',
      'T cell activation'
    )
  ) |>
  dplyr::select(Description, pvalue, Count) |>
  mutate(
    log_pvalue = -log10(pvalue),
    category = "Jurkat"
  )

# Define colors for dot plots
Figure2C_GO_colors <- c(
  "Common" = "#808080",
  "HEK293T" = "#4DBBD5",
  "HepG2" = "#F39B7F",
  "Jurkat" = "#00A087"
)

# Dot plot for Common O-GlcNAc proteins (left side, faceted by Ontology)
Figure2C_common_GO_plot <- ggplot(common_OGlcNAc_GO_selected, aes(x = log_pvalue, y = Description)) +
  geom_segment(
    aes(x = 0, xend = Inf, y = Description, yend = Description),
    linetype = "dashed",
    color = "grey70",
    linewidth = 0.3
  ) +
  geom_point(aes(size = Count), color = "#808080") +
  facet_grid(rows = vars(ONTOLOGY), scales = "free_y", space = "free_y") +
  scale_size_continuous(range = c(1, 4)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Common",
    x = expression(-log[10]('p value')),
    y = NULL,
    size = "Count"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica", color = "black"),
    plot.title = element_text(size = 6, hjust = 0.5),
    axis.text = element_text(color = "black", size = 5),
    axis.text.y = element_text(lineheight = 0.8),
    axis.title = element_text(size = 5),
    strip.background = element_blank(),
    strip.text = element_text(size = 5),
    legend.position = "none"
  )

# Combine unique protein GO terms for faceted plot (right side)
unique_GO_df <- bind_rows(
  unique_HEK293T_GO_selected,
  unique_HepG2_GO_selected,
  unique_Jurkat_GO_selected
) |>
  mutate(
    category = factor(category, levels = c("HEK293T", "HepG2", "Jurkat")),
    Description = str_wrap(Description, width = 30)
  )

# Dot plot for Unique O-GlcNAc proteins (right side, faceted)
Figure2C_unique_GO_plot <- ggplot(unique_GO_df, aes(x = log_pvalue, y = Description)) +
  geom_segment(
    aes(x = 0, xend = Inf, y = Description, yend = Description),
    linetype = "dashed",
    color = "grey70",
    linewidth = 0.3
  ) +
  geom_point(aes(size = Count, color = category)) +
  facet_grid(rows = vars(category), scales = "free_y", space = "free_y") +
  scale_color_manual(values = Figure2C_GO_colors) +
  scale_size_continuous(range = c(1, 4)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Unique",
    x = expression(-log[10]('p value')),
    y = NULL,
    size = "Count"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica", color = "black"),
    plot.title = element_text(size = 6, hjust = 0.5),
    axis.text = element_text(color = "black", size = 5),
    axis.text.y = element_text(lineheight = 0.8),
    axis.title = element_text(size = 5),
    strip.background = element_blank(),
    strip.text = element_text(size = 5),
    legend.position = "none"
  ) +
  guides(color = "none")

# Combine all plots using patchwork
library(patchwork)

# Layout: Common (left) | Donut (center) | Unique (right)
Figure2C_combined <- Figure2C_common_GO_plot + Figure2C_donut + Figure2C_unique_GO_plot +
  plot_layout(widths = c(1, 3, 1))

ggsave(
  filename = paste0(figure_file_path, "Figure2/Figure2C.pdf"),
  plot = Figure2C_combined,
  width = 5, height = 2, units = "in"
)

OGlcNAc_proteins_HepG2 <- OGlyco_HepG2_bonafide |>
  filter(Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859')) |>
  pull(Protein.ID) |>
  unique()

OGlcNAc_proteins_Jurkat <- OGlyco_Jurkat_bonafide |>
  filter(Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859')) |>
  pull(Protein.ID) |>
  unique()

# Extract O-GalNAc proteins from each cell type
OGalNAc_proteins_HEK293T <- OGlyco_HEK293T_bonafide |>
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1) % 326.1339', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968')) |>
  pull(Protein.ID) |>
  unique()

OGalNAc_proteins_HepG2 <- OGlyco_HepG2_bonafide |>
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1) % 326.1339', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968')) |>
  pull(Protein.ID) |>
  unique()

OGalNAc_proteins_Jurkat <- OGlyco_Jurkat_bonafide |>
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1) % 326.1339', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968')) |>
  pull(Protein.ID) |>
  unique()

# Calculate total proteins (union across all 3 cell types)
total_OGlcNAc_proteins <- union(
  union(OGlcNAc_proteins_HEK293T, OGlcNAc_proteins_HepG2),
  OGlcNAc_proteins_Jurkat
)

total_OGalNAc_proteins <- union(
  union(OGalNAc_proteins_HEK293T, OGalNAc_proteins_HepG2),
  OGalNAc_proteins_Jurkat
)

# Calculate common proteins (intersection of all 3 cell types)
common_OGlcNAc_proteins <- intersect(

  intersect(OGlcNAc_proteins_HEK293T, OGlcNAc_proteins_HepG2),
  OGlcNAc_proteins_Jurkat
)

common_OGalNAc_proteins <- intersect(
  intersect(OGalNAc_proteins_HEK293T, OGalNAc_proteins_HepG2),
  OGalNAc_proteins_Jurkat
)

# Calculate counts and percentages
OGlcNAc_total_count <- length(total_OGlcNAc_proteins)
OGlcNAc_common_count <- length(common_OGlcNAc_proteins)
OGlcNAc_common_pct <- OGlcNAc_common_count / OGlcNAc_total_count * 100

OGalNAc_total_count <- length(total_OGalNAc_proteins)
OGalNAc_common_count <- length(common_OGalNAc_proteins)
OGalNAc_common_pct <- OGalNAc_common_count / OGalNAc_total_count * 100

# Create data for radial bar chart
# Background rings (full circle, 0-100%)
# ymin/ymax adjusted to leave hollow center (starting from y=2)
Figure2D_bg_df <- data.frame(
  glycan_type = c("O-GlcNAc", "O-GalNAc"),
  ymin = c(3, 2),
  ymax = c(4, 3),
  xmin = 0,
  xmax = 100,
  alpha_type = "Total"
)

# Filled arcs (common portion) - black color
Figure2D_fill_df <- data.frame(
  glycan_type = c("O-GlcNAc", "O-GalNAc"),
  ymin = c(3, 2),
  ymax = c(4, 3),
  xmin = 0,
  xmax = c(OGlcNAc_common_pct, OGalNAc_common_pct),
  common_count = c(OGlcNAc_common_count, OGalNAc_common_count),
  total_count = c(OGlcNAc_total_count, OGalNAc_total_count),
  alpha_type = "Common"
)

# Labels for total counts (at 12 o'clock position)
Figure2D_total_labels <- data.frame(
  glycan_type = c("O-GlcNAc", "O-GalNAc"),
  x = 0,
  y = c(3.5, 2.5),
  label = c(OGlcNAc_total_count, OGalNAc_total_count)
)

# Labels for common counts (at midpoint of filled arc)
Figure2D_common_labels <- data.frame(
  glycan_type = c("O-GlcNAc", "O-GalNAc"),
  x = c(OGlcNAc_common_pct / 2, OGalNAc_common_pct / 2),
  y = c(3.5, 2.5),
  label = c(OGlcNAc_common_count, OGalNAc_common_count)
)

# Create radial bar chart
Figure2D <- ggplot() +
  # Background rings (full circle) - light glycan colors
  geom_rect(
    data = Figure2D_bg_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = glycan_type, alpha = alpha_type),
    color = 'black', linewidth = 0.2
  ) +
  # Filled arcs (common portion) - darker glycan colors
  geom_rect(
    data = Figure2D_fill_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = glycan_type, alpha = alpha_type),
    color = 'black', linewidth = 0.2
  ) +
  # Total count labels at 12 o'clock
  geom_text(
    data = Figure2D_total_labels,
    aes(x = x, y = y, label = label),
    size = 2,
    color = "black",
    hjust = 1
  ) +
  # Common count labels
  geom_text(
    data = Figure2D_common_labels,
    aes(x = x, y = y, label = label),
    size = 2,
    color = "black"
  ) +
  scale_fill_manual(
    values = colors_glycan,
    labels = c("O-GalNAc", "O-GlcNAc")
  ) +
  scale_alpha_manual(
    values = c("Total" = 0.3, "Common" = 0.9),
    labels = c("Total", "Common")
  ) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 4.5)) +
  coord_polar(theta = "x", start = 0, direction = 1) +
  labs(fill = NULL, alpha = NULL) +
  guides(
    fill = guide_legend(order = 1, ncol = 1),
    alpha = guide_legend(order = 2, ncol = 1)
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.size = unit(0.15, 'cm'),
    legend.key.height = unit(0.15, 'cm'),
    legend.key.width = unit(0.15, 'cm'),
    legend.text = element_text(family = "Helvetica", size = 5),
    legend.spacing.x = unit(0.3, 'cm'),
    legend.spacing.y = unit(0, 'cm'),
    legend.box.spacing = unit(0.05, 'cm'),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 2)
  )

ggsave(
  filename = paste0(figure_file_path, "Figure2/Figure2D.pdf"),
  plot = Figure2D,
  width = 1.5, height = 1.5, units = "in"
)

# Functional enrichment analysis for common O-GlcNAc and O-GalNAc proteins
library(clusterProfiler)
library(org.Hs.eg.db)

# Define universe as union of common O-GlcNAc and O-GalNAc proteins
universe_proteins <- union(common_OGlcNAc_proteins, common_OGalNAc_proteins)

## Gene Ontology analysis
# Common O-GlcNAc proteins GO enrichment
common_OGlcNAc_GO <- enrichGO(
  gene = common_OGlcNAc_proteins,
  OrgDb = org.Hs.eg.db,
  universe = universe_proteins,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  common_OGlcNAc_GO@result,
  paste0(source_file_path, 'enrichment/common_OGlcNAc_GO.csv')
)

# Common O-GalNAc proteins GO enrichment
common_OGalNAc_GO <- enrichGO(
  gene = common_OGalNAc_proteins,
  OrgDb = org.Hs.eg.db,
  universe = universe_proteins,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  common_OGalNAc_GO@result,
  paste0(source_file_path, 'enrichment/common_OGalNAc_GO.csv')
)

## KEGG analysis
# Common O-GlcNAc proteins KEGG enrichment
common_OGlcNAc_KEGG <- enrichKEGG(
  gene = common_OGlcNAc_proteins,
  organism = 'hsa',
  keyType = 'uniprot',
  universe = universe_proteins,
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  common_OGlcNAc_KEGG@result,
  paste0(source_file_path, 'enrichment/common_OGlcNAc_KEGG.csv')
)

# Common O-GalNAc proteins KEGG enrichment
common_OGalNAc_KEGG <- enrichKEGG(
  gene = common_OGalNAc_proteins,
  organism = 'hsa',
  keyType = 'uniprot',
  universe = universe_proteins,
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  common_OGalNAc_KEGG@result,
  paste0(source_file_path, 'enrichment/common_OGalNAc_KEGG.csv')
)

# import functional enrichment results of common O-GlcNAc and O-GalNAc
common_OGlcNAc_GO_selected <- read_csv(
  paste0(source_file_path, 'enrichment/common_OGlcNAc_GO.csv')
) |>
  filter(
    Description %in% c(
      'nucleus',
      'regulation of macromolecule metabolic process',
      'nucleoplasm'
    )
  ) |>
  select(
    Description, pvalue, Count
  ) |>
  mutate(
    log_pvalue = -log10(pvalue),
    glycan_type = "O-GlcNAc"
  )

common_OGalNAc_GO_selected <- read_csv(
  paste0(source_file_path, 'enrichment/common_OGalNAc_GO.csv')
) |>
  filter(
    Description %in% c(
      'endoplasmic reticulum lumen',
      'nuclear outer membrane-endoplasmic reticulum membrane network',
      'protein maturation'
    )
  ) |>
  select(
    Description, pvalue, Count
  ) |>
  mutate(
    log_pvalue = -log10(pvalue),
    glycan_type = "O-GalNAc"
  )

# Combine O-GlcNAc and O-GalNAc selected terms
Figure2D_GO_df <- bind_rows(
  common_OGlcNAc_GO_selected,
  common_OGalNAc_GO_selected
) |>
  mutate(
    glycan_type = factor(glycan_type, levels = c("O-GlcNAc", "O-GalNAc")),
    Description = str_replace(Description, "endoplasmic reticulum", "ER"),
    Description = str_wrap(Description, width = 30)
  )

# Dot plot for functional enrichment
Figure2D_GO_plot <- ggplot(Figure2D_GO_df, aes(x = log_pvalue, y = Description)) +
  # Horizontal dashed lines (from left to right edge)
  geom_segment(
    aes(x = 0, xend = Inf, y = Description, yend = Description),
    linetype = "dashed",
    color = "grey70",
    linewidth = 0.3
  ) +
  # Dots
  geom_point(aes(size = Count, color = glycan_type)) +
  # Facet by glycan type (labels on right, space proportional to content)
  facet_grid(rows = vars(glycan_type), scales = "free_y", space = "free_y") +
  # Colors
  scale_color_manual(values = colors_glycan) +
  # Size scale (smaller circles)
  scale_size_continuous(range = c(1, 4)) +
  # X-axis
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  # Labels
  labs(
    x = expression(-log[10]('p value')),
    y = NULL,
    color = NULL,
    size = "Count"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Helvetica", color = "black"),
    axis.text = element_text(color = "black", size = 6),
    axis.text.y = element_text(lineheight = 0.6),
    axis.title = element_text(size = 6),
    strip.background = element_blank(),
    strip.text = element_text(size = 6),
    legend.key.size = unit(0.2, 'cm'),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 5),
    legend.position = "right"
  ) +
  guides(color = "none")

ggsave(
  filename = paste0(figure_file_path, "Figure2/Figure2D_GO_plot.pdf"),
  plot = Figure2D_GO_plot,
  width = 2.5, height = 1.5, units = "in"
)
