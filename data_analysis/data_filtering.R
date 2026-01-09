# import packages
library(tidyverse)

# data filtering
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


