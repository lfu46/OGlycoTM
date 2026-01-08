# import packages
library(tidyverse)

# source data
source("data_analysis/data_source.R")

# Glycoprotein level quantification ---------------------------------------

# HEK293T
OGlyco_glycoprotein_quant_HEK293T <- OGlyco_HEK293T_bonafide |>
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

write_csv(OGlyco_glycoprotein_quant_HEK293T,
          paste0(source_file_path, "OGlyco_glycoprotein_quant_HEK293T.csv"))

# HepG2
OGlyco_glycoprotein_quant_HepG2 <- OGlyco_HepG2_bonafide |>
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

write_csv(OGlyco_glycoprotein_quant_HepG2,
          paste0(source_file_path, "OGlyco_glycoprotein_quant_HepG2.csv"))

# Jurkat
OGlyco_glycoprotein_quant_Jurkat <- OGlyco_Jurkat_bonafide |>
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

write_csv(OGlyco_glycoprotein_quant_Jurkat,
          paste0(source_file_path, "OGlyco_glycoprotein_quant_Jurkat.csv"))

