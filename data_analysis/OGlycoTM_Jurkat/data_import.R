# import packages
library(tidyverse)

# import OGlcNAc result
OGlcNAc_Jurkat_299_localized <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Jurkat/OGlyco_Search_6/OGlycoTM_Jurkat/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(528\\.2859\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(299\\.1230\\)')
  ) |> 
  filter(!is.na(STYC.299.1230)) |> 
  filter(STYC.299.1230.Best.Localization >= 0.75)

OGlcNAc_Jurkat_299_nolocalized <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Jurkat/OGlyco_Search_6/OGlycoTM_Jurkat/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(528\\.2859\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(299\\.1230\\)')
  ) |> 
  filter(!is.na(STYC.299.1230)) |> 
  filter(STYC.299.1230.Best.Localization < 0.75) |> 
  filter(!str_detect(Peptide, 'C'))

OGlcNAc_Jurkat_528_localized <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Jurkat/OGlyco_Search_6/OGlycoTM_Jurkat/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(528\\.2859\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(299\\.1230\\)')
  ) |> 
  filter(!is.na(STYC.528.2859)) |> 
  filter(STYC.528.2859.Best.Localization >= 0.75)

OGlcNAc_Jurkat_528_nolocalized <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Jurkat/OGlyco_Search_6/OGlycoTM_Jurkat/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(528\\.2859\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(299\\.1230\\)')
  ) |> 
  filter(!is.na(STYC.528.2859)) |> 
  filter(STYC.528.2859.Best.Localization < 0.75) |> 
  filter(!str_detect(Peptide, 'C'))

OGlcNAc_Jurkat_bonafide <- bind_rows(
  OGlcNAc_Jurkat_299_localized,
  OGlcNAc_Jurkat_299_nolocalized,
  OGlcNAc_Jurkat_528_localized,
  OGlcNAc_Jurkat_528_nolocalized
)

write_csv(
  OGlcNAc_Jurkat_bonafide,
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Jurkat/OGlyco_Search_6/OGlcNAc_Jurkat_bonafide.csv'
)

# import OGalNAc result
library(tidyverse)

OGalNAc_Jurkat_326_localized <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Jurkat/OGlyco_Search_7/OGlycoTM_Jurkat/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(555\\.2968\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(326\\.1339\\)')
  ) |> 
  filter(!is.na(STYC.326.1339)) |> 
  filter(STYC.326.1339.Best.Localization >= 0.75)

OGalNAc_Jurkat_326_nolocalized <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Jurkat/OGlyco_Search_7/OGlycoTM_Jurkat/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(555\\.2968\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(326\\.1339\\)')
  ) |> 
  filter(!is.na(STYC.326.1339)) |> 
  filter(STYC.326.1339.Best.Localization < 0.75) |> 
  filter(!str_detect(Peptide, 'C'))

OGalNAc_Jurkat_555_localized <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Jurkat/OGlyco_Search_7/OGlycoTM_Jurkat/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(555\\.2968\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(326\\.1339\\)')
  ) |> 
  filter(!is.na(STYC.555.2968)) |> 
  filter(STYC.555.2968.Best.Localization >= 0.75)

OGalNAc_Jurkat_555_nolocalized <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Jurkat/OGlyco_Search_7/OGlycoTM_Jurkat/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(555\\.2968\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(326\\.1339\\)')
  ) |> 
  filter(!is.na(STYC.555.2968)) |> 
  filter(STYC.555.2968.Best.Localization < 0.75) |> 
  filter(!str_detect(Peptide, 'C'))

OGalNAc_Jurkat_bonafide <- bind_rows(
  OGalNAc_Jurkat_326_localized,
  OGalNAc_Jurkat_326_nolocalized,
  OGalNAc_Jurkat_555_localized,
  OGalNAc_Jurkat_555_nolocalized
)

write_csv(
  OGalNAc_Jurkat_bonafide,
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Jurkat/OGlyco_Search_7/OGalNAc_Jurkat_bonafide.csv'
)
