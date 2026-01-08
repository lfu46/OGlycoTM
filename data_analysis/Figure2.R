# import packages
library(tidyverse)

# Figure 2. Identification of O-GlcNAc and O-GalNAc sites and proteins in HEK293T, HepG2 and Jurkat cells


# Figure 2A ---------------------------------------------------------------
# Indentification of O-GlcNAc and O-GalNAc proteins in three type of cells

# Number of unique O-GlcNAcylated protein
# HEK293T
OGlcNAc_protein_unique_number_HEK293T <- OGlyco_HEK293T_bonafide |> 
  filter(
    Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859')
  ) |> 
  distinct(Protein.ID) |> 
  nrow()

# HepG2
OGlcNAc_protein_unique_number_HepG2 <- OGlyco_HepG2_bonafide |> 
  filter(
    Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859')
  ) |> 
  distinct(Protein.ID) |> 
  nrow()

# Jurkat
OGlcNAc_protein_unique_number_Jurkat <- OGlyco_Jurkat_bonafide |> 
  filter(
    Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859')
  ) |> 
  distinct(Protein.ID) |> 
  nrow()

# Number of unique O-GalNAcylated protein
# HEK293T
OGalNAc_protein_unique_number_HEK293T <- OGlyco_HEK293T_bonafide |> 
  filter(
    Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1) % 326.1339', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968')
  ) |> 
  distinct(Protein.ID) |> 
  nrow()

# HepG2
OGalNAc_protein_unique_number_HepG2 <- OGlyco_HepG2_bonafide |> 
  filter(
    Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1) % 326.1339', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968')
  ) |> 
  distinct(Protein.ID) |> 
  nrow()

# Jurkat
OGalNAc_protein_unique_number_Jurkat <- OGlyco_Jurkat_bonafide |> 
  filter(
    Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1) % 326.1339', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968')
  ) |> 
  distinct(Protein.ID) |> 
  nrow()

# summary of the identification results
Figure2A_df <- data.frame(
  cell_type = rep(c("HEK293T", "HepG2", "Jurkat"), 2),
  glycan_type = factor(
    c(rep("O-GlcNAc", 3), rep("O-GalNAc", 3)),
    levels = c("O-GlcNAc", "O-GalNAc")
  ),
  protein_count = c(
    OGlcNAc_protein_unique_number_HEK293T,
    OGlcNAc_protein_unique_number_HepG2,
    OGlcNAc_protein_unique_number_Jurkat,
    OGalNAc_protein_unique_number_HEK293T,
    OGalNAc_protein_unique_number_HepG2,
    OGalNAc_protein_unique_number_Jurkat
  )
)

# barplot
Figure2A <- ggplot(Figure2A_df, aes(x = cell_type, y = protein_count, fill = glycan_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(
    aes(label = protein_count),
    position = position_dodge(width = 0.9),
    vjust = -0.2,
    size = 2,
    color = "black"
  ) +
  scale_fill_manual(values = colors_glycan) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = NULL,
    y = "No. of unique glycoproteins",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    title = element_text(family = "Helvetica", color = "black", size = 6),
    text = element_text(family = "Helvetica", color = "black"),
    axis.text = element_text(color = "black", size = 6),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6)
  )

ggsave(
  filename = paste0(figure_file_path, "Figure2/Figure2A.pdf"),
  plot = Figure2A,
  width = 2.5, height = 1.5, units = "in"
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
write_csv(OGlyco_site_HEK293T, paste0(source_file_path, "OGlyco_site_HEK293T.csv"))
write_csv(OGlyco_site_HepG2, paste0(source_file_path, "OGlyco_site_HepG2.csv"))
write_csv(OGlyco_site_Jurkat, paste0(source_file_path, "OGlyco_site_Jurkat.csv"))

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

# create list for Venn diagram
OGlcNAc_site_list <- list(
  HEK293T = OGlcNAc_site_HEK293T,
  HepG2 = OGlcNAc_site_HepG2,
  Jurkat = OGlcNAc_site_Jurkat
)

# create euler object (proportional)
OGlcNAc_euler <- euler(OGlcNAc_site_list)

# Venn diagram
Figure2B <- plot(
  OGlcNAc_euler,
  fills = list(fill = colors_cell, alpha = 0.5),
  edges = list(col = "white", lwd = 2),
  labels = list(font = 1, cex = 0.7, nudge = 0.2),
  quantities = list(font = 1, cex = 0.7)
)

# save plot
pdf(paste0(figure_file_path, "Figure2/Figure2B.pdf"), width = 2, height = 1.5)
print(Figure2B)
dev.off()


# Figure 2C ---------------------------------------------------------------
# UniProt, O-GlcNAc Atlas and PhosphoSitePlus annotation of O-GlcNAc sites identified in all cell types
# generate total identified O-GlcNAc proteins
OGlcNAc_protein_total <- bind_rows(
  OGlyco_HEK293T_bonafide |>
    filter(Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859')) |>
    select(Protein.ID),
  OGlyco_HepG2_bonafide |>
    filter(Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859')) |>
    select(Protein.ID),
  OGlyco_Jurkat_bonafide |>
    filter(Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859')) |>
    select(Protein.ID)
) |>
  distinct()

write_csv(
  OGlcNAc_protein_total,
  paste0(source_file_path, 'OGlcNAc_protein_total.csv')
)

# generate total localized O-GlcNAc sites
OGlcNAc_site_total <- bind_rows(
  OGlyco_site_HEK293T |>
    filter(glycan_type == 'O-GlcNAc') |>
    select(site_index),
  OGlyco_site_HepG2 |>
    filter(glycan_type == 'O-GlcNAc') |>
    select(site_index),
  OGlyco_site_Jurkat |>
    filter(glycan_type == 'O-GlcNAc') |>
    select(site_index)
) |>
  distinct()

write_csv(
  OGlcNAc_site_total,
  paste0(source_file_path, 'OGlcNAc_site_total.csv')
)

# process mapping results of UniProt
UniProt_OGlcNAc_mapping <- read_tsv(
  paste0(source_file_path, 'idmapping_2026_01_07_OGlcNAc_protein_total.tsv')
) |>
  select(Entry, Glycosylation) |>
  filter(!is.na(Glycosylation)) |>
  filter(str_detect(Glycosylation, "O-linked \\(GlcNAc\\)"))

# parse glycosylation column to extract O-GlcNAc sites
UniProt_OGlcNAc_sites <- UniProt_OGlcNAc_mapping |>
  mutate(
    # split glycosylation entries by "; CARBOHYD"
    glyco_entries = str_split(Glycosylation, "; (?=CARBOHYD)")
  ) |>
  unnest(glyco_entries) |>
  # keep only O-linked GlcNAc entries (exclude N-linked)
  filter(str_detect(glyco_entries, "O-linked \\(GlcNAc\\)")) |>
  # for alternate sites, require PubMed evidence

  filter(
    !str_detect(glyco_entries, "alternate") |
    str_detect(glyco_entries, "PubMed")
  ) |>
  mutate(
    # extract position number
    position = as.numeric(str_extract(glyco_entries, "(?<=CARBOHYD )\\d+")),
    # extract residue type
    residue_full = str_extract(glyco_entries, "(?<=O-linked \\(GlcNAc\\) )(serine|threonine)"),
    # convert to single letter
    modified_residue = case_when(
      residue_full == "serine" ~ "S",
      residue_full == "threonine" ~ "T"
    ),
    # create site index
    site_index = paste0(Entry, "_", modified_residue, position)
  ) |>
  select(Entry, position, modified_residue, site_index) |>
  distinct()

# process PhosphoSitePlus O-GlcNAc sites
PhosphoSitePlus_OGlcNAc_sites <- read_csv(
  paste0(source_file_path, 'OGlcNAc_site_PhosphoSitePlus.csv'),
  skip = 3
) |>
  filter(ORGANISM == "human") |>
  filter(Ambiguous_Site == 0) |>
  mutate(
    # extract residue (first character, uppercase)
    modified_residue = toupper(str_extract(MOD_RSD, "^[A-Za-z]")),
    # extract position number
    position = as.numeric(str_extract(MOD_RSD, "\\d+")),
    # create site index
    site_index = paste0(ACC_ID, "_", modified_residue, position)
  ) |>
  select(ACC_ID, position, modified_residue, site_index) |>
  distinct()

# process OGlcNAc Atlas sites
OGlcNAcAtlas_sites <- read_csv(
  paste0(source_file_path, 'OGlcNAcAtlas_unambiguous_sites_20251208.csv')
) |>
  filter(species == "human") |>
  mutate(
    site_index = paste0(accession, "_", site_residue, position_in_protein)
  ) |>
  select(accession, position_in_protein, site_residue, site_index) |>
  distinct()

# combine all database annotations
database_annotated_sites <- bind_rows(
  UniProt_OGlcNAc_sites |> select(site_index),
  PhosphoSitePlus_OGlcNAc_sites |> select(site_index),
  OGlcNAcAtlas_sites |> select(site_index)
) |>
  distinct()

# calculate annotated vs not annotated sites
annotated_count <- OGlcNAc_site_total |>
  filter(site_index %in% database_annotated_sites$site_index) |>
  nrow()

not_annotated_count <- OGlcNAc_site_total |>
  filter(!site_index %in% database_annotated_sites$site_index) |>
  nrow()

# create dataframe for stacked bar plot
Figure2C_df <- data.frame(
  x_label = c("O-GlcNAc"),
  annotation_status = factor(
    c("Annotated", "Not annotated"),
    levels = c("Not annotated", "Annotated")
  ),
  site_count = c(annotated_count, not_annotated_count)
)

# stacked bar plot
Figure2C <- ggplot(Figure2C_df, aes(x = x_label, y = site_count, fill = annotation_status)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = NULL,
    y = "No. of unique O-GlcNAc sites",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    title = element_text(family = "Helvetica", color = "black", size = 6),
    text = element_text(family = "Helvetica", color = "black"),
    axis.text = element_text(color = "black", size = 6),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6)
  )

ggsave(
  filename = paste0(figure_file_path, "Figure2/Figure2C.pdf"),
  plot = Figure2C,
  width = 1.8, height = 1.5, units = "in"
)


# Figure 2D -----------------------------------------------------
# Radial bar chart showing percentage of common O-GlcNAc and O-GalNAc proteins
# Functional enrichment analysis of common O-GlcNAc and O-GalNAc proteins
# Two concentric rings: outer = O-GlcNAc, inner = O-GalNAc

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
  paste0(source_file_path, 'common_OGlcNAc_GO.csv')
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
  paste0(source_file_path, 'common_OGalNAc_GO.csv')
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
  paste0(source_file_path, 'common_OGlcNAc_KEGG.csv')
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
  paste0(source_file_path, 'common_OGalNAc_KEGG.csv')
)

# import functional enrichment results of common O-GlcNAc and O-GalNAc
common_OGlcNAc_GO_selected <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/common_OGlcNAc_GO.csv'
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
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/common_OGalNAc_GO.csv'
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


# Figure 2E ---------------------------------------------------------------
# Radial bar chart showing percentage of unique O-GlcNAc and O-GalNAc proteins in HEK293T
# Functional enrichment analysis for unique O-GlcNAc and O-GalNAc proteins in HEK293T

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

# Calculate unique HEK293T proteins (only in HEK293T, not in HepG2 or Jurkat)
unique_OGlcNAc_HEK293T <- setdiff(
  setdiff(OGlcNAc_proteins_HEK293T, OGlcNAc_proteins_HepG2),
  OGlcNAc_proteins_Jurkat
)

unique_OGalNAc_HEK293T <- setdiff(
  setdiff(OGalNAc_proteins_HEK293T, OGalNAc_proteins_HepG2),
  OGalNAc_proteins_Jurkat
)

# Calculate counts and percentages
OGlcNAc_total_count <- length(total_OGlcNAc_proteins)
OGlcNAc_unique_HEK293T_count <- length(unique_OGlcNAc_HEK293T)
OGlcNAc_unique_HEK293T_pct <- OGlcNAc_unique_HEK293T_count / OGlcNAc_total_count * 100

OGalNAc_total_count <- length(total_OGalNAc_proteins)
OGalNAc_unique_HEK293T_count <- length(unique_OGalNAc_HEK293T)
OGalNAc_unique_HEK293T_pct <- OGalNAc_unique_HEK293T_count / OGalNAc_total_count * 100

# Create data for radial bar chart
# Background rings (full circle, 0-100%)
Figure2E_bg_df <- data.frame(
  glycan_type = c("O-GlcNAc", "O-GalNAc"),
  ymin = c(3, 2),
  ymax = c(4, 3),
  xmin = 0,
  xmax = 100,
  alpha_type = "Total"
)

# Filled arcs (unique HEK293T portion)
Figure2E_fill_df <- data.frame(
  glycan_type = c("O-GlcNAc", "O-GalNAc"),
  ymin = c(3, 2),
  ymax = c(4, 3),
  xmin = 0,
  xmax = c(OGlcNAc_unique_HEK293T_pct, OGalNAc_unique_HEK293T_pct),
  unique_count = c(OGlcNAc_unique_HEK293T_count, OGalNAc_unique_HEK293T_count),
  total_count = c(OGlcNAc_total_count, OGalNAc_total_count),
  alpha_type = "Unique"
)

# Labels for total counts (at 12 o'clock position)
Figure2E_total_labels <- data.frame(
  glycan_type = c("O-GlcNAc", "O-GalNAc"),
  x = 0,
  y = c(3.5, 2.5),
  label = c(OGlcNAc_total_count, OGalNAc_total_count)
)

# Labels for unique counts (at midpoint of filled arc)
Figure2E_unique_labels <- data.frame(
  glycan_type = c("O-GlcNAc", "O-GalNAc"),
  x = c(OGlcNAc_unique_HEK293T_pct / 2, OGalNAc_unique_HEK293T_pct / 2),
  y = c(3.5, 2.5),
  label = c(OGlcNAc_unique_HEK293T_count, OGalNAc_unique_HEK293T_count)
)

# Create radial bar chart
Figure2E <- ggplot() +
  # Background rings (full circle) - light grey with black border
  geom_rect(
    data = Figure2E_bg_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "grey90",
    color = "black", linewidth = 0.2
  ) +
  # Filled arcs (unique portion) - glycan colors without border
  geom_rect(
    data = Figure2E_fill_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = glycan_type),
    color = "black", linewidth = 0.2
  ) +
  # Total count labels at 12 o'clock
  geom_text(
    data = Figure2E_total_labels,
    aes(x = x, y = y, label = label),
    size = 2,
    color = "black",
    hjust = 1
  ) +
  # Unique count labels
  geom_text(
    data = Figure2E_unique_labels,
    aes(x = x, y = y, label = label),
    size = 2,
    color = "black"
  ) +
  scale_fill_manual(
    values = colors_glycan,
    labels = c("O-GalNAc", "O-GlcNAc")
  ) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 4.5)) +
  coord_polar(theta = "x", start = 0, direction = 1) +
  labs(fill = NULL) +
  guides(
    fill = guide_legend(nrow = 1)
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.size = unit(0.15, 'cm'),
    legend.key.height = unit(0.15, 'cm'),
    legend.key.width = unit(0.15, 'cm'),
    legend.text = element_text(family = "Helvetica", size = 5),
    legend.spacing.x = unit(0.1, 'cm'),
    legend.spacing.y = unit(0, 'cm'),
    legend.box.spacing = unit(0.05, 'cm'),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 2)
  )

ggsave(
  filename = paste0(figure_file_path, "Figure2/Figure2E.pdf"),
  plot = Figure2E,
  width = 1.5, height = 1.5, units = "in"
)

# Functional enrichment analysis for unique O-GlcNAc and O-GalNAc proteins in HEK293T
library(clusterProfiler)
library(org.Hs.eg.db)

## Gene Ontology analysis for unique O-GlcNAc proteins in HEK293T
unique_OGlcNAc_HEK293T_GO <- enrichGO(
  gene = unique_OGlcNAc_HEK293T,
  OrgDb = org.Hs.eg.db,
  universe = total_OGlcNAc_proteins,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  unique_OGlcNAc_HEK293T_GO@result,
  paste0(source_file_path, 'unique_OGlcNAc_HEK293T_GO.csv')
)

## Gene Ontology analysis for unique O-GalNAc proteins in HEK293T
unique_OGalNAc_HEK293T_GO <- enrichGO(
  gene = unique_OGalNAc_HEK293T,
  OrgDb = org.Hs.eg.db,
  universe = total_OGalNAc_proteins,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  unique_OGalNAc_HEK293T_GO@result,
  paste0(source_file_path, 'unique_OGalNAc_HEK293T_GO.csv')
)

## KEGG analysis for unique O-GlcNAc proteins in HEK293T
unique_OGlcNAc_HEK293T_KEGG <- enrichKEGG(
  gene = unique_OGlcNAc_HEK293T,
  organism = 'hsa',
  keyType = 'uniprot',
  universe = total_OGlcNAc_proteins,
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  unique_OGlcNAc_HEK293T_KEGG@result,
  paste0(source_file_path, 'unique_OGlcNAc_HEK293T_KEGG.csv')
)

## KEGG analysis for unique O-GalNAc proteins in HEK293T
unique_OGalNAc_HEK293T_KEGG <- enrichKEGG(
  gene = unique_OGalNAc_HEK293T,
  organism = 'hsa',
  keyType = 'uniprot',
  universe = total_OGalNAc_proteins,
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  unique_OGalNAc_HEK293T_KEGG@result,
  paste0(source_file_path, 'unique_OGalNAc_HEK293T_KEGG.csv')
)

# import functional enrichment results of unique O-GlcNAc and O-GalNAc in HEK293T
unique_OGlcNAc_HEK293T_GO_selected <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/unique_OGlcNAc_HEK293T_GO.csv'
) |>
  filter(
    Description %in% c(
      'chaperone-mediated protein folding',
      'cell cycle process',
      'epithelial cell proliferation'
    )
  ) |>
  select(
    Description, pvalue, Count
  ) |>
  mutate(
    log_pvalue = -log10(pvalue),
    glycan_type = "O-GlcNAc"
  )

unique_OGalNAc_HEK293T_GO_selected <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/unique_OGalNAc_HEK293T_GO.csv'
) |>
  filter(
    Description %in% c(
      'cytosol',
      'membraneless organelle'
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
Figure2E_GO_df <- bind_rows(
  unique_OGlcNAc_HEK293T_GO_selected,
  unique_OGalNAc_HEK293T_GO_selected
) |>
  mutate(
    glycan_type = factor(glycan_type, levels = c("O-GlcNAc", "O-GalNAc")),
    Description = str_wrap(Description, width = 30)
  )

# Dot plot for functional enrichment
Figure2E_GO_plot <- ggplot(Figure2E_GO_df, aes(x = log_pvalue, y = Description)) +
  # Horizontal dashed lines
  geom_segment(
    aes(x = 0, xend = Inf, y = Description, yend = Description),
    linetype = "dashed",
    color = "grey70",
    linewidth = 0.3
  ) +
  # Dots
  geom_point(aes(size = Count, color = glycan_type)) +
  # Facet by glycan type (space proportional to content)
  facet_grid(rows = vars(glycan_type), scales = "free_y", space = "free_y") +
  # Colors
  scale_color_manual(values = colors_glycan) +
  # Size scale
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
  filename = paste0(figure_file_path, "Figure2/Figure2E_GO_plot.pdf"),
  plot = Figure2E_GO_plot,
  width = 2.5, height = 1.5, units = "in"
)


# Figure 2F ---------------------------------------------------------------
# Functional enrichment analysis for unique O-GlcNAc and O-GalNAc proteins in each cell type
library(clusterProfiler)
library(org.Hs.eg.db)

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

# Extract unique O-GlcNAc proteins for each cell type
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

# Extract unique O-GalNAc proteins for each cell type
unique_OGalNAc_HEK293T <- setdiff(
  setdiff(OGalNAc_proteins_HEK293T, OGalNAc_proteins_HepG2),
  OGalNAc_proteins_Jurkat
)

unique_OGalNAc_HepG2 <- setdiff(
  setdiff(OGalNAc_proteins_HepG2, OGalNAc_proteins_HEK293T),
  OGalNAc_proteins_Jurkat
)

unique_OGalNAc_Jurkat <- setdiff(
  setdiff(OGalNAc_proteins_Jurkat, OGalNAc_proteins_HEK293T),
  OGalNAc_proteins_HepG2
)

# Define universe as total O-GlcNAc and O-GalNAc proteins
total_OGlcNAc_proteins <- union(
  union(OGlcNAc_proteins_HEK293T, OGlcNAc_proteins_HepG2),
  OGlcNAc_proteins_Jurkat
)

total_OGalNAc_proteins <- union(
  union(OGalNAc_proteins_HEK293T, OGalNAc_proteins_HepG2),
  OGalNAc_proteins_Jurkat
)

## Gene Ontology analysis for unique O-GlcNAc proteins
# Unique O-GlcNAc HEK293T GO enrichment
unique_OGlcNAc_HEK293T_GO <- enrichGO(
  gene = unique_OGlcNAc_HEK293T,
  OrgDb = org.Hs.eg.db,
  universe = total_OGlcNAc_proteins,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  unique_OGlcNAc_HEK293T_GO@result,
  paste0(source_file_path, 'unique_OGlcNAc_HEK293T_GO.csv')
)

# Unique O-GlcNAc HepG2 GO enrichment
unique_OGlcNAc_HepG2_GO <- enrichGO(
  gene = unique_OGlcNAc_HepG2,
  OrgDb = org.Hs.eg.db,
  universe = total_OGlcNAc_proteins,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  unique_OGlcNAc_HepG2_GO@result,
  paste0(source_file_path, 'unique_OGlcNAc_HepG2_GO.csv')
)

# Unique O-GlcNAc Jurkat GO enrichment
unique_OGlcNAc_Jurkat_GO <- enrichGO(
  gene = unique_OGlcNAc_Jurkat,
  OrgDb = org.Hs.eg.db,
  universe = total_OGlcNAc_proteins,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  unique_OGlcNAc_Jurkat_GO@result,
  paste0(source_file_path, 'unique_OGlcNAc_Jurkat_GO.csv')
)

## Gene Ontology analysis for unique O-GalNAc proteins
# Unique O-GalNAc HEK293T GO enrichment
unique_OGalNAc_HEK293T_GO <- enrichGO(
  gene = unique_OGalNAc_HEK293T,
  OrgDb = org.Hs.eg.db,
  universe = total_OGalNAc_proteins,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  unique_OGalNAc_HEK293T_GO@result,
  paste0(source_file_path, 'unique_OGalNAc_HEK293T_GO.csv')
)

# Unique O-GalNAc HepG2 GO enrichment
unique_OGalNAc_HepG2_GO <- enrichGO(
  gene = unique_OGalNAc_HepG2,
  OrgDb = org.Hs.eg.db,
  universe = total_OGalNAc_proteins,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  unique_OGalNAc_HepG2_GO@result,
  paste0(source_file_path, 'unique_OGalNAc_HepG2_GO.csv')
)

# Unique O-GalNAc Jurkat GO enrichment
unique_OGalNAc_Jurkat_GO <- enrichGO(
  gene = unique_OGalNAc_Jurkat,
  OrgDb = org.Hs.eg.db,
  universe = total_OGalNAc_proteins,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  unique_OGalNAc_Jurkat_GO@result,
  paste0(source_file_path, 'unique_OGalNAc_Jurkat_GO.csv')
)

## KEGG analysis for unique O-GlcNAc proteins
# Unique O-GlcNAc HEK293T KEGG enrichment
unique_OGlcNAc_HEK293T_KEGG <- enrichKEGG(
  gene = unique_OGlcNAc_HEK293T,
  organism = 'hsa',
  keyType = 'uniprot',
  universe = total_OGlcNAc_proteins,
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  unique_OGlcNAc_HEK293T_KEGG@result,
  paste0(source_file_path, 'unique_OGlcNAc_HEK293T_KEGG.csv')
)

# Unique O-GlcNAc HepG2 KEGG enrichment
unique_OGlcNAc_HepG2_KEGG <- enrichKEGG(
  gene = unique_OGlcNAc_HepG2,
  organism = 'hsa',
  keyType = 'uniprot',
  universe = total_OGlcNAc_proteins,
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  unique_OGlcNAc_HepG2_KEGG@result,
  paste0(source_file_path, 'unique_OGlcNAc_HepG2_KEGG.csv')
)

# Unique O-GlcNAc Jurkat KEGG enrichment
unique_OGlcNAc_Jurkat_KEGG <- enrichKEGG(
  gene = unique_OGlcNAc_Jurkat,
  organism = 'hsa',
  keyType = 'uniprot',
  universe = total_OGlcNAc_proteins,
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  unique_OGlcNAc_Jurkat_KEGG@result,
  paste0(source_file_path, 'unique_OGlcNAc_Jurkat_KEGG.csv')
)

## KEGG analysis for unique O-GalNAc proteins
# Unique O-GalNAc HEK293T KEGG enrichment
unique_OGalNAc_HEK293T_KEGG <- enrichKEGG(
  gene = unique_OGalNAc_HEK293T,
  organism = 'hsa',
  keyType = 'uniprot',
  universe = total_OGalNAc_proteins,
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  unique_OGalNAc_HEK293T_KEGG@result,
  paste0(source_file_path, 'unique_OGalNAc_HEK293T_KEGG.csv')
)

# Unique O-GalNAc HepG2 KEGG enrichment
unique_OGalNAc_HepG2_KEGG <- enrichKEGG(
  gene = unique_OGalNAc_HepG2,
  organism = 'hsa',
  keyType = 'uniprot',
  universe = total_OGalNAc_proteins,
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  unique_OGalNAc_HepG2_KEGG@result,
  paste0(source_file_path, 'unique_OGalNAc_HepG2_KEGG.csv')
)

# Unique O-GalNAc Jurkat KEGG enrichment
unique_OGalNAc_Jurkat_KEGG <- enrichKEGG(
  gene = unique_OGalNAc_Jurkat,
  organism = 'hsa',
  keyType = 'uniprot',
  universe = total_OGalNAc_proteins,
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  unique_OGalNAc_Jurkat_KEGG@result,
  paste0(source_file_path, 'unique_OGalNAc_Jurkat_KEGG.csv')
)

