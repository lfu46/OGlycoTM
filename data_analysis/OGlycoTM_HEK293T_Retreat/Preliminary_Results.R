# import packages
library(tidyverse)

# import HEK293T Preliminary Results
Preliminary_Results <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OGlyco_HEK293T_Preliminary_Results/OGlycoTM_EThcD_OPair_TMT_Search_1/OGlyco_HEK293T_Preliminary_Results/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  filter(!is.na(O.Pair.Score)) |> 
  select(
    Spectrum,
    Peptide, Peptide.Length, Charge, Observed.Mass, Observed.M.Z, Delta.Mass, Hyperscore, Nextscore,
    Number.of.Enzymatic.Termini, Number.of.Missed.Cleavages, Protein.Start, Protein.End, 
    Assigned.Modifications, O.Pair.Score:Paired.Scan.Num, Parent.Scan.Number,
    Protein.ID:Protein.Description
  )

# data filtering
# exclude glycopeptides with localized sites on cysteine resiudes
# exclude glycopeptides that contains cysteine in their sequences but lack localized sites
Preliminary_Results_localized <- Preliminary_Results |> 
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

Preliminary_Results_nolocalized <- Preliminary_Results |> 
  filter(!Confidence.Level %in% c('Level1', 'Level1b')) |> 
  filter(!str_detect(Peptide, 'C'))

Preliminary_Results_bonafide <- bind_rows(
  Preliminary_Results_localized,
  Preliminary_Results_nolocalized
)

write_csv(
  Preliminary_Results_bonafide,
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OGlyco_HEK293T_Preliminary_Results/Preliminary_Results_bonafide.csv'
)

# import filtered results
library(tidyverse)

Preliminary_Results_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OGlyco_HEK293T_Preliminary_Results/Preliminary_Results_bonafide.csv'
)

# check number of total glycoPSM and unique glycoprotein
# O-GlcNAc
Total_GlycoPSM_OGlcNAc <- Preliminary_Results_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)', 'HexNAt(1)TMT6plex(1)')) |> 
  nrow()

Unique_Glycoprotein_OGlcNAc <- Preliminary_Results_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)', 'HexNAt(1)TMT6plex(1)')) |> 
  distinct(Protein.ID) |> 
  nrow()

GlycoPSM_Glycoprotein_Results_OGlcNAc <- tibble(
  category = c('Total GlycoPSM', 'Unique Glycoprotein'),
  number = c(5651, 566)
)

# O-GalNAc
Total_GlycoPSM_OGalNAc <- Preliminary_Results_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1)', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1)')) |> 
  nrow()

Unique_Glycoprotein_OGalNAc <- Preliminary_Results_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1)', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1)')) |> 
  distinct(Protein.ID) |> 
  nrow()

GlycoPSM_Glycoprotein_Results_OGalNAc <- tibble(
  category = c('Total GlycoPSM', 'Unique Glycoprotein'),
  number = c(550, 129)
)

# check number of localized and non-localized glycoPSM
# O-GlcNAc
Localized_GlycoPSM_OGlcNAc <- Preliminary_Results_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)', 'HexNAt(1)TMT6plex(1)')) |> 
  filter(Confidence.Level %in% c('Level1', 'Level1b')) |> 
  nrow()

Nonlocalized_GlycoPSM_OGlcNAc <- Preliminary_Results_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)', 'HexNAt(1)TMT6plex(1)')) |> 
  filter(!Confidence.Level %in% c('Level1', 'Level1b')) |> 
  nrow()

Localized_Nonlocalized_Results_OGlcNAc <- tibble(
  category = c('Localized', 'Nonlocalized'),
  number = c(552, 5099)
)

# O-GalNAc
Localized_GlycoPSM_OGalNAc <- Preliminary_Results_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1)', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1)')) |> 
  filter(Confidence.Level %in% c('Level1', 'Level1b')) |> 
  nrow()

Nonlocalized_GlycoPSM_OGalNAc <- Preliminary_Results_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1)', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1)')) |> 
  filter(!Confidence.Level %in% c('Level1', 'Level1b')) |> 
  nrow()

Localized_Nonlocalized_Results_OGalNAc <- tibble(
  category = c('Localized', 'Nonlocalized'),
  number = c(137, 413)
)



