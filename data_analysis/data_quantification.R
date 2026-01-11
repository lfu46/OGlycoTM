# import packages
library(tidyverse)

# source data
source("data_source.R")

# ==============================================================================
# O-GlcNAc Filtered Datasets
# ==============================================================================

# Filter bonafide data for O-GlcNAc modifications only
OGlcNAc_HEK293T <- OGlyco_HEK293T_bonafide |>
  filter(Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859'))

OGlcNAc_HepG2 <- OGlyco_HepG2_bonafide |>
  filter(Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859'))

OGlcNAc_Jurkat <- OGlyco_Jurkat_bonafide |>
  filter(Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859'))

# Save O-GlcNAc filtered datasets
write_csv(OGlcNAc_HEK293T, paste0(source_file_path, "filtered/OGlcNAc_HEK293T.csv"))
write_csv(OGlcNAc_HepG2, paste0(source_file_path, "filtered/OGlcNAc_HepG2.csv"))
write_csv(OGlcNAc_Jurkat, paste0(source_file_path, "filtered/OGlcNAc_Jurkat.csv"))


# ==============================================================================
# O-GlcNAc Protein Level Quantification
# ==============================================================================

# HEK293T
OGlcNAc_protein_quant_HEK293T <- OGlcNAc_HEK293T |>
  group_by(Protein.ID) |>
  summarize(
    Entry.Name = first(Entry.Name),
    Gene = first(Gene),
    Protein.Description = first(Protein.Description),
    Intensity.Tuni_1 = sum(Intensity.Tuni_1),
    Intensity.Tuni_2 = sum(Intensity.Tuni_2),
    Intensity.Tuni_3 = sum(Intensity.Tuni_3),
    Intensity.Ctrl_4 = sum(Intensity.Ctrl_4),
    Intensity.Ctrl_5 = sum(Intensity.Ctrl_5),
    Intensity.Ctrl_6 = sum(Intensity.Ctrl_6)
  )

write_csv(OGlcNAc_protein_quant_HEK293T,
          paste0(source_file_path, "quantification/OGlcNAc_protein_quant_HEK293T.csv"))

# HepG2
OGlcNAc_protein_quant_HepG2 <- OGlcNAc_HepG2 |>
  group_by(Protein.ID) |>
  summarize(
    Entry.Name = first(Entry.Name),
    Gene = first(Gene),
    Protein.Description = first(Protein.Description),
    Intensity.Tuni_1 = sum(Intensity.Tuni_1),
    Intensity.Tuni_2 = sum(Intensity.Tuni_2),
    Intensity.Tuni_3 = sum(Intensity.Tuni_3),
    Intensity.Ctrl_4 = sum(Intensity.Ctrl_4),
    Intensity.Ctrl_5 = sum(Intensity.Ctrl_5),
    Intensity.Ctrl_6 = sum(Intensity.Ctrl_6)
  )

write_csv(OGlcNAc_protein_quant_HepG2,
          paste0(source_file_path, "quantification/OGlcNAc_protein_quant_HepG2.csv"))

# Jurkat
OGlcNAc_protein_quant_Jurkat <- OGlcNAc_Jurkat |>
  group_by(Protein.ID) |>
  summarize(
    Entry.Name = first(Entry.Name),
    Gene = first(Gene),
    Protein.Description = first(Protein.Description),
    Intensity.Tuni_1 = sum(Intensity.Tuni_1),
    Intensity.Tuni_2 = sum(Intensity.Tuni_2),
    Intensity.Tuni_3 = sum(Intensity.Tuni_3),
    Intensity.Ctrl_4 = sum(Intensity.Ctrl_4),
    Intensity.Ctrl_5 = sum(Intensity.Ctrl_5),
    Intensity.Ctrl_6 = sum(Intensity.Ctrl_6)
  )

write_csv(OGlcNAc_protein_quant_Jurkat,
          paste0(source_file_path, "quantification/OGlcNAc_protein_quant_Jurkat.csv"))


# ==============================================================================
# O-GlcNAc Site Filtered Datasets
# ==============================================================================

# Filter site data for O-GlcNAc modifications only
OGlcNAc_site_HEK293T <- OGlyco_site_HEK293T |>
  filter(glycan_type == "O-GlcNAc")

OGlcNAc_site_HepG2 <- OGlyco_site_HepG2 |>
  filter(glycan_type == "O-GlcNAc")

OGlcNAc_site_Jurkat <- OGlyco_site_Jurkat |>
  filter(glycan_type == "O-GlcNAc")

# Save O-GlcNAc site filtered datasets
write_csv(OGlcNAc_site_HEK293T, paste0(source_file_path, "site/OGlcNAc_site_HEK293T.csv"))
write_csv(OGlcNAc_site_HepG2, paste0(source_file_path, "site/OGlcNAc_site_HepG2.csv"))
write_csv(OGlcNAc_site_Jurkat, paste0(source_file_path, "site/OGlcNAc_site_Jurkat.csv"))


# ==============================================================================
# O-GlcNAc Site Level Quantification
# ==============================================================================

# HEK293T
OGlcNAc_site_quant_HEK293T <- OGlcNAc_site_HEK293T |>
  group_by(site_index) |>
  summarize(
    Protein.ID = first(Protein.ID),
    Entry.Name = first(Entry.Name),
    Gene = first(Gene),
    Protein.Description = first(Protein.Description),
    modified_residue = first(modified_residue),
    site_number = first(site_number),
    glycan_type = first(glycan_type),
    Intensity.Tuni_1 = sum(Intensity.Tuni_1),
    Intensity.Tuni_2 = sum(Intensity.Tuni_2),
    Intensity.Tuni_3 = sum(Intensity.Tuni_3),
    Intensity.Ctrl_4 = sum(Intensity.Ctrl_4),
    Intensity.Ctrl_5 = sum(Intensity.Ctrl_5),
    Intensity.Ctrl_6 = sum(Intensity.Ctrl_6)
  )

write_csv(OGlcNAc_site_quant_HEK293T,
          paste0(source_file_path, "quantification/OGlcNAc_site_quant_HEK293T.csv"))

# HepG2
OGlcNAc_site_quant_HepG2 <- OGlcNAc_site_HepG2 |>
  group_by(site_index) |>
  summarize(
    Protein.ID = first(Protein.ID),
    Entry.Name = first(Entry.Name),
    Gene = first(Gene),
    Protein.Description = first(Protein.Description),
    modified_residue = first(modified_residue),
    site_number = first(site_number),
    glycan_type = first(glycan_type),
    Intensity.Tuni_1 = sum(Intensity.Tuni_1),
    Intensity.Tuni_2 = sum(Intensity.Tuni_2),
    Intensity.Tuni_3 = sum(Intensity.Tuni_3),
    Intensity.Ctrl_4 = sum(Intensity.Ctrl_4),
    Intensity.Ctrl_5 = sum(Intensity.Ctrl_5),
    Intensity.Ctrl_6 = sum(Intensity.Ctrl_6)
  )

write_csv(OGlcNAc_site_quant_HepG2,
          paste0(source_file_path, "quantification/OGlcNAc_site_quant_HepG2.csv"))

# Jurkat
OGlcNAc_site_quant_Jurkat <- OGlcNAc_site_Jurkat |>
  group_by(site_index) |>
  summarize(
    Protein.ID = first(Protein.ID),
    Entry.Name = first(Entry.Name),
    Gene = first(Gene),
    Protein.Description = first(Protein.Description),
    modified_residue = first(modified_residue),
    site_number = first(site_number),
    glycan_type = first(glycan_type),
    Intensity.Tuni_1 = sum(Intensity.Tuni_1),
    Intensity.Tuni_2 = sum(Intensity.Tuni_2),
    Intensity.Tuni_3 = sum(Intensity.Tuni_3),
    Intensity.Ctrl_4 = sum(Intensity.Ctrl_4),
    Intensity.Ctrl_5 = sum(Intensity.Ctrl_5),
    Intensity.Ctrl_6 = sum(Intensity.Ctrl_6)
  )

write_csv(OGlcNAc_site_quant_Jurkat,
          paste0(source_file_path, "quantification/OGlcNAc_site_quant_Jurkat.csv"))


# ==============================================================================
# Whole Proteome Protein Level Quantification
# ==============================================================================

# Read filtered whole proteome data
WP_HEK293T_filtered <- read_csv(
  paste0(source_file_path, 'filtered/WP_HEK293T_filtered.csv')
)
WP_HepG2_filtered <- read_csv(
  paste0(source_file_path, 'filtered/WP_HepG2_filtered.csv')
)
WP_Jurkat_filtered <- read_csv(
  paste0(source_file_path, 'filtered/WP_Jurkat_filtered.csv')
)

# HEK293T
WP_protein_quant_HEK293T <- WP_HEK293T_filtered |>
  group_by(UniProt_Accession) |>
  summarize(
    Gene.Symbol = first(Gene.Symbol),
    Protein.ID = first(Protein.ID),
    Protein.MWT.kDa. = first(Protein.MWT.kDa.),
    Annotation = first(Annotation),
    Intensity.Tuni_1 = sum(Sn.126),
    Intensity.Tuni_2 = sum(Sn.127),
    Intensity.Tuni_3 = sum(Sn.128),
    Intensity.Ctrl_4 = sum(Sn.129),
    Intensity.Ctrl_5 = sum(Sn.130),
    Intensity.Ctrl_6 = sum(Sn.131)
  )

write_csv(WP_protein_quant_HEK293T,
          paste0(source_file_path, "quantification/WP_protein_quant_HEK293T.csv"))

# HepG2
WP_protein_quant_HepG2 <- WP_HepG2_filtered |>
  group_by(UniProt_Accession) |>
  summarize(
    Gene.Symbol = first(Gene.Symbol),
    Protein.ID = first(Protein.ID),
    Protein.MWT.kDa. = first(Protein.MWT.kDa.),
    Annotation = first(Annotation),
    Intensity.Tuni_1 = sum(Sn.126),
    Intensity.Tuni_2 = sum(Sn.127),
    Intensity.Tuni_3 = sum(Sn.128),
    Intensity.Ctrl_4 = sum(Sn.129),
    Intensity.Ctrl_5 = sum(Sn.130),
    Intensity.Ctrl_6 = sum(Sn.131)
  )

write_csv(WP_protein_quant_HepG2,
          paste0(source_file_path, "quantification/WP_protein_quant_HepG2.csv"))

# Jurkat
WP_protein_quant_Jurkat <- WP_Jurkat_filtered |>
  group_by(UniProt_Accession) |>
  summarize(
    Gene.Symbol = first(Gene.Symbol),
    Protein.ID = first(Protein.ID),
    Protein.MWT.kDa. = first(Protein.MWT.kDa.),
    Annotation = first(Annotation),
    Intensity.Tuni_1 = sum(Sn.126),
    Intensity.Tuni_2 = sum(Sn.127),
    Intensity.Tuni_3 = sum(Sn.128),
    Intensity.Ctrl_4 = sum(Sn.129),
    Intensity.Ctrl_5 = sum(Sn.130),
    Intensity.Ctrl_6 = sum(Sn.131)
  )

write_csv(WP_protein_quant_Jurkat,
          paste0(source_file_path, "quantification/WP_protein_quant_Jurkat.csv"))

