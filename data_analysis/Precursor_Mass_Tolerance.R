# import packages
library(tidyverse)

# import HEK293T Preliminary Results
Preliminary_Results_Processed <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OGlyco_HEK293T_Preliminary_Results/OGlycoTM_EThcD_OPair_TMT_Search_1/OGlyco_HEK293T_Preliminary_Results/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  mutate(
    PPM = Delta.Mass/Calculated.Peptide.Mass
  )

write_csv(
  Preliminary_Results_Processed,
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OGlyco_HEK293T_Preliminary_Results/Preliminary_Results_Processed.csv'
)

