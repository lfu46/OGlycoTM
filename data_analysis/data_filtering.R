# import packages
library(tidyverse)

# ==============================================================================
# OGlyco Data Filtering
# ==============================================================================

# HEK293T
# exclude glycopeptides with localized sites on cysteine resiudes
OGlyco_HEK293T_localized <- OGlyco_HEK293T_raw |> 
  filter(Confidence.Level %in% c('Level1', 'Level1b')) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(555\\.2968\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(528\\.2859\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(326\\.1339\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(299\\.1230\\)')
  )

# exclude glycopeptides that contains cysteine in their sequences but lack localized sites
OGlyco_HEK293T_nolocalized <- OGlyco_HEK293T_raw |> 
  filter(!Confidence.Level %in% c('Level1', 'Level1b')) |> 
  filter(!str_detect(Peptide, 'C'))

OGlyco_HEK293T_bonafide <- bind_rows(
  OGlyco_HEK293T_localized,
  OGlyco_HEK293T_nolocalized
)

write_csv(
  OGlyco_HEK293T_bonafide,
  file = paste0(source_file_path, 'filtered/OGlyco_HEK293T_bonafide.csv')
)

# HepG2
# exclude glycopeptides with localized sites on cysteine resiudes
OGlyco_HepG2_localized <- OGlyco_HepG2_raw |> 
  filter(Confidence.Level %in% c('Level1', 'Level1b')) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(555\\.2968\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(528\\.2859\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(326\\.1339\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(299\\.1230\\)')
  )

# exclude glycopeptides that contains cysteine in their sequences but lack localized sites
OGlyco_HepG2_nolocalized <- OGlyco_HepG2_raw |> 
  filter(!Confidence.Level %in% c('Level1', 'Level1b')) |> 
  filter(!str_detect(Peptide, 'C'))

OGlyco_HepG2_bonafide <- bind_rows(
  OGlyco_HepG2_localized,
  OGlyco_HepG2_nolocalized
)

write_csv(
  OGlyco_HepG2_bonafide,
  file = paste0(source_file_path, 'filtered/OGlyco_HepG2_bonafide.csv')
)

# Jurkat
# exclude glycopeptides with localized sites on cysteine resiudes
OGlyco_Jurkat_localized <- OGlyco_Jurkat_raw |> 
  filter(Confidence.Level %in% c('Level1', 'Level1b')) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(555\\.2968\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(528\\.2859\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(326\\.1339\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(299\\.1230\\)')
  )

# exclude glycopeptides that contains cysteine in their sequences but lack localized sites
OGlyco_Jurkat_nolocalized <- OGlyco_Jurkat_raw |> 
  filter(!Confidence.Level %in% c('Level1', 'Level1b')) |> 
  filter(!str_detect(Peptide, 'C'))

OGlyco_Jurkat_bonafide <- bind_rows(
  OGlyco_Jurkat_localized,
  OGlyco_Jurkat_nolocalized
)

write_csv(
  OGlyco_Jurkat_bonafide,
  file = paste0(source_file_path, 'filtered/OGlyco_Jurkat_bonafide.csv')
)

# ==============================================================================
# Whole Proteome Data Filtering
# ==============================================================================

# Read raw whole proteome data
WP_HEK293T_raw <- read_csv(
  paste0(source_file_path, 'raw/WP_HEK293T_raw.csv')
)
WP_HepG2_raw <- read_csv(
  paste0(source_file_path, 'raw/WP_HepG2_raw.csv')
)
WP_Jurkat_raw <- read_csv(
  paste0(source_file_path, 'raw/WP_Jurkat_raw.csv')
)

# HEK293T
# Filter: remove decoys (##), remove contaminants, XCorr > 1.2, PPM -10 to 10, all Sn > 5
# Extract UniProt_Accession from Reference column
WP_HEK293T_filtered <- WP_HEK293T_raw |>
  filter(!str_detect(Reference, '^##')) |>
  filter(!str_detect(Reference, regex('contaminant', ignore_case = TRUE))) |>
  filter(!str_detect(Protein.ID, '^##')) |>
  filter(!str_detect(Protein.ID, regex('contaminant', ignore_case = TRUE))) |>
  filter(XCorr > 1.2) |>
  filter(PPM > -10, PPM < 10) |>
  filter(Sn.126 > 5, Sn.127 > 5, Sn.128 > 5, Sn.129 > 5, Sn.130 > 5, Sn.131 > 5) |>
  mutate(
    UniProt_Accession = str_extract(Reference, '(?<=\\|)[A-Z0-9]+(?=\\|)')
  )

write_csv(
  WP_HEK293T_filtered,
  file = paste0(source_file_path, 'filtered/WP_HEK293T_filtered.csv')
)

# HepG2
WP_HepG2_filtered <- WP_HepG2_raw |>
  filter(!str_detect(Reference, '^##')) |>
  filter(!str_detect(Reference, regex('contaminant', ignore_case = TRUE))) |>
  filter(!str_detect(Protein.ID, '^##')) |>
  filter(!str_detect(Protein.ID, regex('contaminant', ignore_case = TRUE))) |>
  filter(XCorr > 1.2) |>
  filter(PPM > -10, PPM < 10) |>
  filter(Sn.126 > 5, Sn.127 > 5, Sn.128 > 5, Sn.129 > 5, Sn.130 > 5, Sn.131 > 5) |>
  mutate(
    UniProt_Accession = str_extract(Reference, '(?<=\\|)[A-Z0-9]+(?=\\|)')
  )

write_csv(
  WP_HepG2_filtered,
  file = paste0(source_file_path, 'filtered/WP_HepG2_filtered.csv')
)

# Jurkat
WP_Jurkat_filtered <- WP_Jurkat_raw |>
  filter(!str_detect(Reference, '^##')) |>
  filter(!str_detect(Reference, regex('contaminant', ignore_case = TRUE))) |>
  filter(!str_detect(Protein.ID, '^##')) |>
  filter(!str_detect(Protein.ID, regex('contaminant', ignore_case = TRUE))) |>
  filter(XCorr > 1.2) |>
  filter(PPM > -10, PPM < 10) |>
  filter(Sn.126 > 5, Sn.127 > 5, Sn.128 > 5, Sn.129 > 5, Sn.130 > 5, Sn.131 > 5) |>
  mutate(
    UniProt_Accession = str_extract(Reference, '(?<=\\|)[A-Z0-9]+(?=\\|)')
  )

write_csv(
  WP_Jurkat_filtered,
  file = paste0(source_file_path, 'filtered/WP_Jurkat_filtered.csv')
)

