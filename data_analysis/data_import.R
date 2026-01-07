# import packages
library(tidyverse)

# import OGlyco results
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
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/OGlyco_HEK293T_raw.csv'
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
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/OGlyco_HepG2_raw.csv'
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
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/OGlyco_Jurkat_raw.csv'
)
