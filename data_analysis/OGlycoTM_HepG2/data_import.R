## import packages
library(tidyverse)

# import results from psm.tsv file
OGlyco_HepG2_raw <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HepG2/EThcD_OPair_Search_1/OGlycoTM_HepG2/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  filter(!is.na(Total.Glycan.Composition)) |> 
  select(
    Spectrum, 
    Peptide, Peptide.Length, Charge, 
    Observed.Mass, Observed.M.Z, Delta.Mass, Hyperscore, Nextscore,
    Number.of.Enzymatic.Termini, Number.of.Missed.Cleavages, Protein.Start, Protein.End, 
    Assigned.Modifications, O.Pair.Score:Site.Probabilities,
    Ratio.138.144 = ..138.144.Ratio, 
    Has.N.Glyc.Sequon, Paired.Scan.Num, Parent.Scan.Number,
    Protein.ID:Protein.Description
  )

## data filtering
## exclude glycopeptides with localized sites on cysteine residues
## exclude glycopeptides that contain cysteine in their sequences but lack localized sites
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

OGlyco_HepG2_nolocalized <- OGlyco_HepG2_raw |> 
  filter(!Confidence.Level %in% c('Level1', 'Level1b')) |> 
  filter(!str_detect(Peptide, 'C'))

OGlyco_HepG2_bonafide <- bind_rows(
  OGlyco_HepG2_localized,
  OGlyco_HepG2_nolocalized
)

write_csv(
  OGlyco_HepG2_bonafide,
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HepG2/EThcD_OPair_Search_1/OGlyco_HepG2_bonafide.csv'
)

# import filtered results
library(tidyverse)

OGlyco_HepG2_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HepG2/EThcD_OPair_Search_1/OGlyco_HepG2_bonafide.csv'
)

# check preliminary results
# O-GlcNAc
Total_GlycoPSM_OGlcNAc <- OGlyco_HepG2_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)', 'HexNAt(1)TMT6plex(1)')) |> 
  nrow()

Unique_Glycoprotein_OGlcNAc <- OGlyco_HepG2_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)', 'HexNAt(1)TMT6plex(1)')) |> 
  distinct(Protein.ID) |> 
  nrow()

GlycoPSM_Glycoprotein_Results_OGlcNAc <- tibble(
  category = c('Total GlycoPSM', 'Unique Glycoprotein'),
  number = c(4963, 594)
)

# O-GalNAc
Total_GlycoPSM_OGalNAc <- OGlyco_HepG2_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1)', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1)')) |> 
  nrow()

Unique_Glycoprotein_OGalNAc <- OGlyco_HepG2_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1)', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1)')) |> 
  distinct(Protein.ID) |> 
  nrow()

GlycoPSM_Glycoprotein_Results_OGalNAc <- tibble(
  category = c('Total GlycoPSM', 'Unique Glycoprotein'),
  number = c(429, 128)
)

# check the total glycoPSM group by raw file
# O-GlcNAc
Total_GlycoPSM_OGlcNAc_by_raw_file <- OGlyco_HepG2_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)', 'HexNAt(1)TMT6plex(1)')) |> 
  mutate(
    raw_file = sub("\\..*", "", Spectrum)
  ) |> 
  group_by(raw_file) |> 
  count() |> 
  ungroup()

# O-GalNAc
Total_GlycoPSM_OGalNAc_by_raw_file <- OGlyco_HepG2_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1)', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1)')) |> 
  mutate(
    raw_file = sub("\\..*", "", Spectrum)
  ) |> 
  group_by(raw_file) |> 
  count() |> 
  ungroup()
