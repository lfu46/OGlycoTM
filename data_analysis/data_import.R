# import packages
library(tidyverse)

# ==============================================================================
# OGlyco Data Import
# ==============================================================================

# HEK293T
OGlyco_HEK293T_raw <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_HEK293T/OGlyco/EThcD_OPair_TMT_Search/OGlyco/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  filter(Is.Decoy == FALSE & Is.Contaminant == FALSE) |> 
  filter(!is.na(Total.Glycan.Composition)) |> 
  select(
    Spectrum,
    Peptide, Peptide.Length, Charge, 
    Observed.Mass, Calculated.Peptide.Mass, Observed.M.Z, Delta.Mass, Hyperscore, Nextscore,
    Number.of.Enzymatic.Termini, Number.of.Missed.Cleavages, Protein.Start, Protein.End, 
    Assigned.Modifications, O.Pair.Score:Total.Glycan.Composition, Confidence.Level:Site.Probabilities,
    Ratio.138.144 = ..138.144.Ratio,
    Has.N.Glyc.Sequon, Paired.Scan.Num,
    Protein.ID:Protein.Description,
    Intensity.Tuni_1:Intensity.Ctrl_6
  )

write_csv(
  OGlyco_HEK293T_raw,
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/raw/OGlyco_HEK293T_raw.csv'
)

# HepG2
OGlyco_HepG2_raw <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_HepG2/OGlyco/EThcD_OPair_TMT_Search/OGlyco/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  filter(Is.Decoy == FALSE & Is.Contaminant == FALSE) |> 
  filter(!is.na(Total.Glycan.Composition)) |> 
  select(
    Spectrum,
    Peptide, Peptide.Length, Charge, 
    Observed.Mass, Calculated.Peptide.Mass, Observed.M.Z, Delta.Mass, Hyperscore, Nextscore,
    Number.of.Enzymatic.Termini, Number.of.Missed.Cleavages, Protein.Start, Protein.End, 
    Assigned.Modifications, O.Pair.Score:Total.Glycan.Composition, Confidence.Level:Site.Probabilities,
    Ratio.138.144 = ..138.144.Ratio,
    Has.N.Glyc.Sequon, Paired.Scan.Num,
    Protein.ID:Protein.Description,
    Intensity.Tuni_1:Intensity.Ctrl_6
  )

write_csv(
  OGlyco_HepG2_raw,
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/raw/OGlyco_HepG2_raw.csv'
)

# Jurkat
OGlyco_Jurkat_raw <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_Jurkat/OGlyco/EThcD_OPair_TMT_Search/OGlyco/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  filter(Is.Decoy == FALSE & Is.Contaminant == FALSE) |> 
  filter(!is.na(Total.Glycan.Composition)) |> 
  select(
    Spectrum,
    Peptide, Peptide.Length, Charge, 
    Observed.Mass, Calculated.Peptide.Mass, Observed.M.Z, Delta.Mass, Hyperscore, Nextscore,
    Number.of.Enzymatic.Termini, Number.of.Missed.Cleavages, Protein.Start, Protein.End, 
    Assigned.Modifications, O.Pair.Score:Total.Glycan.Composition, Confidence.Level:Site.Probabilities,
    Ratio.138.144 = ..138.144.Ratio,
    Has.N.Glyc.Sequon, Paired.Scan.Num,
    Protein.ID:Protein.Description,
    Intensity.Tuni_1:Intensity.Ctrl_6
  )

write_csv(
  OGlyco_Jurkat_raw,
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/raw/OGlyco_Jurkat_raw.csv'
)

# ==============================================================================
# Whole Proteome Data Import
# ==============================================================================

# HEK293T
WP_HEK293T_raw <- read_csv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/raw/ronghuwulab_1768133089_HEK293T.csv',
  col_names = TRUE,
  name_repair = 'universal'
) |>
  rename(
    Sn.126 = `..126.Sn`,
    Sn.127 = `..127.Sn`,
    Sn.128 = `..128.Sn`,
    Sn.129 = `..129.Sn`,
    Sn.130 = `..130.Sn`,
    Sn.131 = `..131.Sn`
  )

write_csv(
  WP_HEK293T_raw,
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/raw/WP_HEK293T_raw.csv'
)

# HepG2
WP_HepG2_raw <- read_csv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/raw/ronghuwulab_1768132256_HepG2.csv',
  col_names = TRUE,
  name_repair = 'universal'
) |>
  rename(
    Sn.126 = `..126.Sn`,
    Sn.127 = `..127.Sn`,
    Sn.128 = `..128.Sn`,
    Sn.129 = `..129.Sn`,
    Sn.130 = `..130.Sn`,
    Sn.131 = `..131.Sn`
  )

write_csv(
  WP_HepG2_raw,
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/raw/WP_HepG2_raw.csv'
)

# Jurkat
WP_Jurkat_raw <- read_csv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/raw/ronghuwulab_1768132802_Jurkat.csv',
  col_names = TRUE,
  name_repair = 'universal'
) |>
  rename(
    Sn.126 = `..126.Sn`,
    Sn.127 = `..127.Sn`,
    Sn.128 = `..128.Sn`,
    Sn.129 = `..129.Sn`,
    Sn.130 = `..130.Sn`,
    Sn.131 = `..131.Sn`
  )

write_csv(
  WP_Jurkat_raw,
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/raw/WP_Jurkat_raw.csv'
)
