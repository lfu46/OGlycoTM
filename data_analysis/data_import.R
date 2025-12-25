# import packages
library(tidyverse)

# import HEK293T OGlyco results
OGlyco_HEK293T <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OGlycoTM_HEK293T/OGlyco/OGlyco_OPair_Search/OGlyco/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  filter(!is.na(O.Pair.Score)) |> 
  select(
    Spectrum,
    Peptide, Peptide.Length, Charge, 
    Observed.Mass, Calculated.Peptide.Mass, Observed.M.Z, Delta.Mass, Hyperscore, Nextscore,
    Number.of.Enzymatic.Termini, Number.of.Missed.Cleavages, Protein.Start, Protein.End, 
    Assigned.Modifications, O.Pair.Score:Site.Probabilities, 
    Ratio.138.144 = ..138.144.Ratio,
    Has.N.Glyc.Sequon, Paired.Scan.Num,
    Protein.ID:Protein.Description
  )

# data filtering
# exclude glycopeptides with localized sites on cysteine resiudes
# exclude glycopeptides that contains cysteine in their sequences but lack localized sites
OGlyco_HEK293T_localized <- OGlyco_HEK293T |> 
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

OGlyco_HEK293T_nolocalized <- OGlyco_HEK293T |> 
  filter(!Confidence.Level %in% c('Level1', 'Level1b')) |> 
  filter(!str_detect(Peptide, 'C'))

OGlyco_HEK293T_bonafide <- bind_rows(
  OGlyco_HEK293T_localized,
  OGlyco_HEK293T_nolocalized
)

write_csv(
  OGlyco_HEK293T_bonafide,
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/data_source/OGlyco_HEK293T_bonafide.csv'
)

# import bonafide glyco results
library(tidyverse)

OGlyco_HEK293T_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/data_source/OGlyco_HEK293T_bonafide.csv'
)

# define the mapping for Delta.Mass ranges to glycan compositions
create_glycan_composition <- function(delta_mass_int) {
  if (delta_mass_int %in% 555:557) {
    return("HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1)")
  } else if (delta_mass_int %in% 528:530) {
    return("HexNAt(1)TMT6plex(1)")
  } else if (delta_mass_int %in% 326:328) {
    return("HexNAt(1)GAO_Methoxylamine(1)")
  } else if (delta_mass_int %in% 299:301) {
    return("HexNAt(1)")
  } else {
    return(NA)
  }
}

# define the mapping for glycan compositions to delta mass corrections
get_mass_correction <- function(glycan_comp) {
  switch(glycan_comp,
         "HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1)" = 555.2968,
         "HexNAt(1)TMT6plex(1)" = 528.2859,
         "HexNAt(1)GAO_Methoxylamine(1)" = 326.1339,
         "HexNAt(1)" = 299.1230,
         NA
  )
}

# reassign glycan composition
OGlyco_HEK293T_bonafide_reassigned <- OGlyco_HEK293T_bonafide |> 
  mutate(
    Total.Glycan.Composition = case_when(
      O.Pair.Score == "No paired scan" ~ sapply(as.integer(Delta.Mass), create_glycan_composition),
      TRUE ~ Total.Glycan.Composition
    ),
    Delta.Mass.Reassigned = case_when(
      O.Pair.Score == "No paired scan" ~ Delta.Mass - sapply(Total.Glycan.Composition, get_mass_correction),
      TRUE ~ Delta.Mass
    )
  ) |> 
  mutate(
    PPM = if_else(
      !is.na(Calculated.Peptide.Mass) & Calculated.Peptide.Mass > 0,
      Delta.Mass.Reassigned / Calculated.Peptide.Mass,
      NA_real_
    )
  ) |>
  select(
    Spectrum,
    Peptide, Peptide.Length, Charge, 
    Observed.Mass, Calculated.Peptide.Mass, Observed.M.Z, Delta.Mass, Delta.Mass.Reassigned, PPM, Hyperscore, Nextscore,
    Number.of.Enzymatic.Termini:Protein.Description
  )

write_csv(
  OGlyco_HEK293T_bonafide_reassigned,
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/data_source/OGlyco_HEK293T_bonafide_reassigned.csv'
)

