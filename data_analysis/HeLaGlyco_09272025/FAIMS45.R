## import packages
library(tidyverse)

# FAIMS45
FAIMS45_Rep1 <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS45/OGlycoTM_EThcD_Search/FAIMS45_1/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  filter(!is.na(O.Pair.Score)) |> 
  select(
    Peptide, Peptide.Length, Charge, Observed.Mass, Observed.M.Z, Delta.Mass, Hyperscore, Nextscore,
    Number.of.Enzymatic.Termini, Number.of.Missed.Cleavages, Protein.Start, Protein.End, 
    Assigned.Modifications, O.Pair.Score:Paired.Scan.Num, 
    Protein.ID:Protein.Description
  )

FAIMS45_Rep2 <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS45/OGlycoTM_EThcD_Search/FAIMS45_2/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  filter(!is.na(O.Pair.Score)) |> 
  select(
    Peptide, Peptide.Length, Charge, Observed.Mass, Observed.M.Z, Delta.Mass, Hyperscore, Nextscore,
    Number.of.Enzymatic.Termini, Number.of.Missed.Cleavages, Protein.Start, Protein.End, 
    Assigned.Modifications, O.Pair.Score:Paired.Scan.Num, 
    Protein.ID:Protein.Description
  )

FAIMS45_Rep3 <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS45/OGlycoTM_EThcD_Search/FAIMS45_3/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  filter(!is.na(O.Pair.Score)) |> 
  select(
    Peptide, Peptide.Length, Charge, Observed.Mass, Observed.M.Z, Delta.Mass, Hyperscore, Nextscore,
    Number.of.Enzymatic.Termini, Number.of.Missed.Cleavages, Protein.Start, Protein.End, 
    Assigned.Modifications, O.Pair.Score:Paired.Scan.Num, 
    Protein.ID:Protein.Description
  )

## data filtering
## exclude glycopeptides with localized sites on cysteine residues
## exclude glycopeptides that contain cysteine in their sequences but lack localized sites
# FAIMS45
# Rep1
FAIMS45_Rep1_localized <- FAIMS45_Rep1 |> 
  filter(Confidence.Level %in% c('Level1', 'Level1b')) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(326\\.1339\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(299\\.1230\\)')
  )

FAIMS45_Rep1_nolocalized <- FAIMS45_Rep1 |> 
  filter(!Confidence.Level %in% c('Level1', 'Level1b')) |> 
  filter(!str_detect(Peptide, 'C'))

FAIMS45_Rep1_bonafide <- bind_rows(
  FAIMS45_Rep1_localized,
  FAIMS45_Rep1_nolocalized
)

write_csv(
  FAIMS45_Rep1_bonafide,
  file = '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS45/OGlycoTM_EThcD_Search/FAIMS45_Rep1_bonafide.csv'
)

# Rep2
FAIMS45_Rep2_localized <- FAIMS45_Rep2 |> 
  filter(Confidence.Level %in% c('Level1', 'Level1b')) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(326\\.1339\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(299\\.1230\\)')
  )

FAIMS45_Rep2_nolocalized <- FAIMS45_Rep2 |> 
  filter(!Confidence.Level %in% c('Level1', 'Level1b')) |> 
  filter(!str_detect(Peptide, 'C'))

FAIMS45_Rep2_bonafide <- bind_rows(
  FAIMS45_Rep2_localized,
  FAIMS45_Rep2_nolocalized
)

write_csv(
  FAIMS45_Rep2_bonafide,
  file = '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS45/OGlycoTM_EThcD_Search/FAIMS45_Rep2_bonafide.csv'
)

# Rep3
FAIMS45_Rep3_localized <- FAIMS45_Rep3 |> 
  filter(Confidence.Level %in% c('Level1', 'Level1b')) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(326\\.1339\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(299\\.1230\\)')
  )

FAIMS45_Rep3_nolocalized <- FAIMS45_Rep3 |> 
  filter(!Confidence.Level %in% c('Level1', 'Level1b')) |> 
  filter(!str_detect(Peptide, 'C'))

FAIMS45_Rep3_bonafide <- bind_rows(
  FAIMS45_Rep3_localized,
  FAIMS45_Rep3_nolocalized
)

write_csv(
  FAIMS45_Rep3_bonafide,
  file = '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS45/OGlycoTM_EThcD_Search/FAIMS45_Rep3_bonafide.csv'
)
