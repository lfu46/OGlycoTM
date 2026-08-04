# data_quantification_OGalNAc.R
# O-GalNAc Protein Level Quantification
# Following the same pipeline as O-GlcNAc in data_quantification.R

library(tidyverse)

# source data
source("data_source.R")

# ==============================================================================
# O-GalNAc Filtered Datasets
# ==============================================================================

# Filter bonafide data for O-GalNAc modifications only
OGalNAc_HEK293T <- OGlyco_HEK293T_bonafide |>
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1) % 326.1339',
                                          'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968'))

OGalNAc_HepG2 <- OGlyco_HepG2_bonafide |>
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1) % 326.1339',
                                          'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968'))

OGalNAc_Jurkat <- OGlyco_Jurkat_bonafide |>
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1) % 326.1339',
                                          'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968'))

# Save O-GalNAc filtered datasets
write_csv(OGalNAc_HEK293T, paste0(source_file_path, "filtered/OGalNAc_HEK293T.csv"))
write_csv(OGalNAc_HepG2, paste0(source_file_path, "filtered/OGalNAc_HepG2.csv"))
write_csv(OGalNAc_Jurkat, paste0(source_file_path, "filtered/OGalNAc_Jurkat.csv"))

cat("O-GalNAc PSMs filtered:\n")
cat("  HEK293T:", nrow(OGalNAc_HEK293T), "PSMs,", n_distinct(OGalNAc_HEK293T$Protein.ID), "proteins\n")
cat("  HepG2:", nrow(OGalNAc_HepG2), "PSMs,", n_distinct(OGalNAc_HepG2$Protein.ID), "proteins\n")
cat("  Jurkat:", nrow(OGalNAc_Jurkat), "PSMs,", n_distinct(OGalNAc_Jurkat$Protein.ID), "proteins\n")

# ==============================================================================
# O-GalNAc Protein Level Quantification
# ==============================================================================

# HEK293T
OGalNAc_protein_quant_HEK293T <- OGalNAc_HEK293T |>
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

write_csv(OGalNAc_protein_quant_HEK293T,
          paste0(source_file_path, "quantification/OGalNAc_protein_quant_HEK293T.csv"))

# HepG2
OGalNAc_protein_quant_HepG2 <- OGalNAc_HepG2 |>
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

write_csv(OGalNAc_protein_quant_HepG2,
          paste0(source_file_path, "quantification/OGalNAc_protein_quant_HepG2.csv"))

# Jurkat
OGalNAc_protein_quant_Jurkat <- OGalNAc_Jurkat |>
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

write_csv(OGalNAc_protein_quant_Jurkat,
          paste0(source_file_path, "quantification/OGalNAc_protein_quant_Jurkat.csv"))

cat("\nO-GalNAc protein quantification:\n")
cat("  HEK293T:", nrow(OGalNAc_protein_quant_HEK293T), "proteins\n")
cat("  HepG2:", nrow(OGalNAc_protein_quant_HepG2), "proteins\n")
cat("  Jurkat:", nrow(OGalNAc_protein_quant_Jurkat), "proteins\n")
