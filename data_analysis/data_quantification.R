# import packages
library(tidyverse)

# source data
source("data_source.R")

# O-GlcNAc filtered datasets ----------------------------------------------

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


# O-GlcNAc protein level quantification -----------------------------------

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


# O-GlcNAc site filtered datasets -----------------------------------------

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


# O-GlcNAc site level quantification --------------------------------------

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

