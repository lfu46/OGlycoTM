# import packages
library(tidyverse)

# import HEK293T OGlyco results
OGlyco_EThcD_HEK293T <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OGlycoTM_HEK293T/OGlyco/OGlyco_EThcD_Search/OGlyco/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  filter(!is.na(Total.Glycan.Composition)) |> 
  mutate(
    PPM = Delta.Mass/Calculated.Peptide.Mass
  ) |>
  select(
    Spectrum,
    Peptide, Peptide.Length, Charge, 
    Observed.Mass, Calculated.Peptide.Mass, Observed.M.Z, Delta.Mass, PPM, Hyperscore, Nextscore,
    Number.of.Enzymatic.Termini, Number.of.Missed.Cleavages, Protein.Start, Protein.End, 
    Assigned.Modifications, Total.Glycan.Composition:Glycan.q.value,
    Protein.ID:Protein.Description
  )
  
OGlyco_EThcD_HEK293T_bonafide <- OGlyco_EThcD_HEK293T |> 
  filter(!str_detect(Peptide, 'C'))

write_csv(
  OGlyco_EThcD_HEK293T_bonafide,
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/data_source/OGlyco_EThcD_HEK293T_bonafide.csv'
)

# import filtered results
library(tidyverse)

OGlyco_EThcD_HEK293T_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/data_source/OGlyco_EThcD_HEK293T_bonafide.csv'
)

# check number of total glycoPSM and unique glycoprotein
# O-GlcNAc
Total_GlycoPSM_OGlcNAc <- OGlyco_EThcD_HEK293T_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859')) |> 
  nrow()

Unique_Glycoprotein_OGlcNAc <- OGlyco_EThcD_HEK293T_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859')) |> 
  distinct(Protein.ID) |> 
  nrow()

GlycoPSM_Glycoprotein_Results_OGlcNAc <- tibble(
  category = c('Total GlycoPSM', 'Unique Glycoprotein'),
  number = c(9059, 699)
)

# O-GalNAc
Total_GlycoPSM_OGalNAc <- OGlyco_EThcD_HEK293T_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968', 'HexNAt(1)GAO_Methoxylamine(1) % 326.1339')) |> 
  nrow()

Unique_Glycoprotein_OGalNAc <- OGlyco_EThcD_HEK293T_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968', 'HexNAt(1)GAO_Methoxylamine(1) % 326.1339')) |> 
  distinct(Protein.ID) |> 
  nrow()

GlycoPSM_Glycoprotein_Results_OGalNAc <- tibble(
  category = c('Total GlycoPSM', 'Unique Glycoprotein'),
  number = c(878, 139)
)


