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


# Figure 2D ---------------------------------------------------------------
# Percentage of the common and unique glycoprotein in each cell type

# Define colors for overlap categories
colors_overlap <- c("Common" = "#4DBBD5", "Shared" = "grey", "Unique" = "#F39B7F")

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

# Function to calculate overlap categories for each cell type
calculate_overlap_categories <- function(proteins_target, proteins_other1, proteins_other2) {
  # Common: in all 3 cell types
  common <- intersect(intersect(proteins_target, proteins_other1), proteins_other2)

  # Unique: only in target cell type
  unique_proteins <- setdiff(setdiff(proteins_target, proteins_other1), proteins_other2)

  # Shared: in exactly 2 cell types (target + one other)
  in_target_and_other1 <- intersect(proteins_target, proteins_other1)
  in_target_and_other2 <- intersect(proteins_target, proteins_other2)
  shared <- setdiff(union(in_target_and_other1, in_target_and_other2), common)

  total <- length(proteins_target)

  data.frame(
    category = c("Common", "Shared", "Unique"),
    count = c(length(common), length(shared), length(unique_proteins)),
    percentage = c(length(common), length(shared), length(unique_proteins)) / total * 100
  )
}

# Calculate for O-GlcNAc
OGlcNAc_HEK293T_overlap <- calculate_overlap_categories(
  OGlcNAc_proteins_HEK293T, OGlcNAc_proteins_HepG2, OGlcNAc_proteins_Jurkat
) |> mutate(cell_type = "HEK293T", glycan_type = "O-GlcNAc")

OGlcNAc_HepG2_overlap <- calculate_overlap_categories(
  OGlcNAc_proteins_HepG2, OGlcNAc_proteins_HEK293T, OGlcNAc_proteins_Jurkat
) |> mutate(cell_type = "HepG2", glycan_type = "O-GlcNAc")

OGlcNAc_Jurkat_overlap <- calculate_overlap_categories(
  OGlcNAc_proteins_Jurkat, OGlcNAc_proteins_HEK293T, OGlcNAc_proteins_HepG2
) |> mutate(cell_type = "Jurkat", glycan_type = "O-GlcNAc")

# Calculate for O-GalNAc
OGalNAc_HEK293T_overlap <- calculate_overlap_categories(
  OGalNAc_proteins_HEK293T, OGalNAc_proteins_HepG2, OGalNAc_proteins_Jurkat
) |> mutate(cell_type = "HEK293T", glycan_type = "O-GalNAc")

OGalNAc_HepG2_overlap <- calculate_overlap_categories(
  OGalNAc_proteins_HepG2, OGalNAc_proteins_HEK293T, OGalNAc_proteins_Jurkat
) |> mutate(cell_type = "HepG2", glycan_type = "O-GalNAc")

OGalNAc_Jurkat_overlap <- calculate_overlap_categories(
  OGalNAc_proteins_Jurkat, OGalNAc_proteins_HEK293T, OGalNAc_proteins_HepG2
) |> mutate(cell_type = "Jurkat", glycan_type = "O-GalNAc")

# Combine all data
Figure2D_df <- bind_rows(
  OGlcNAc_HEK293T_overlap,
  OGlcNAc_HepG2_overlap,
  OGlcNAc_Jurkat_overlap,
  OGalNAc_HEK293T_overlap,
  OGalNAc_HepG2_overlap,
  OGalNAc_Jurkat_overlap
) |>
  mutate(
    category = factor(category, levels = c("Unique", "Shared", "Common")),
    glycan_type = factor(glycan_type, levels = c("O-GlcNAc", "O-GalNAc")),
    cell_type = factor(cell_type, levels = c("HEK293T", "HepG2", "Jurkat"))
  )

# Stacked bar plot
Figure2D <- ggplot(Figure2D_df, aes(x = glycan_type, y = percentage, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = ifelse(category == "Shared", "", count)),
    position = position_stack(vjust = 0.5),
    size = 2,
    color = "black"
  ) +
  facet_wrap(~cell_type, nrow = 1) +
  scale_fill_manual(values = colors_overlap, breaks = c("Common", "Unique")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = NULL,
    y = "Percentage (%)",
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
    legend.title = element_text(size = 6),
    strip.background = element_blank(),
    strip.text = element_text(size = 6)
  )

ggsave(
  filename = paste0(figure_file_path, "Figure2/Figure2D.pdf"),
  plot = Figure2D,
  width = 2.5, height = 1.5, units = "in"
)


# Figure 2E ---------------------------------------------------------------


