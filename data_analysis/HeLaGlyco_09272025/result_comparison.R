# import packags
library(tidyverse)

## import bona fide OGlyco results
# FAIMS45
FAIMS45_Rep1_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS45/OGlycoTM_EThcD_Search/FAIMS45_Rep1_bonafide.csv'
)
FAIMS45_Rep2_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS45/OGlycoTM_EThcD_Search/FAIMS45_Rep2_bonafide.csv'
)
FAIMS45_Rep3_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS45/OGlycoTM_EThcD_Search/FAIMS45_Rep3_bonafide.csv'
)

# FAIMS50
FAIMS50_Rep1_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS50/OGlycoTM_EThcD_Search/FAIMS50_Rep1_bonafide.csv'
)
FAIMS50_Rep2_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS50/OGlycoTM_EThcD_Search/FAIMS50_Rep2_bonafide.csv'
)
FAIMS50_Rep3_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS50/OGlycoTM_EThcD_Search/FAIMS50_Rep3_bonafide.csv'
)

# FAIMS4045
FAIMS4045_Rep1_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS4045/OGlycoTM_EThcD_Search/FAIMS4045_Rep1_bonafide.csv'
)
FAIMS4045_Rep2_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS4045/OGlycoTM_EThcD_Search/FAIMS4045_Rep2_bonafide.csv'
)
FAIMS4045_Rep3_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS4045/OGlycoTM_EThcD_Search/FAIMS4045_Rep3_bonafide.csv'
)

# FAIMS4050
FAIMS4050_Rep1_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS4050/OGlycoTM_EThcD_Search/FAIMS4050_Rep1_bonafide.csv'
)
FAIMS4050_Rep2_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS4050/OGlycoTM_EThcD_Search/FAIMS4050_Rep2_bonafide.csv'
)
FAIMS4050_Rep3_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS4050/OGlycoTM_EThcD_Search/FAIMS4050_Rep3_bonafide.csv'
)

# FAIMS404550
FAIMS404550_Rep1_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS404550/OGlycoTM_EThcD_Search/FAIMS404550_Rep1_bonafide.csv'
)
FAIMS404550_Rep2_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS404550/OGlycoTM_EThcD_Search/FAIMS404550_Rep2_bonafide.csv'
)
FAIMS404550_Rep3_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS404550/OGlycoTM_EThcD_Search/FAIMS404550_Rep3_bonafide.csv'
)

# FAIMS455065
FAIMS455065_Rep1_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS455065/OGlycoTM_EThcD_Search/FAIMS455065_Rep1_bonafide.csv'
)
FAIMS455065_Rep2_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS455065/OGlycoTM_EThcD_Search/FAIMS455065_Rep2_bonafide.csv'
)
FAIMS455065_Rep3_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/FAIMS455065/OGlycoTM_EThcD_Search/FAIMS455065_Rep3_bonafide.csv'
)

# noFAIMS_OT
noFAIMS_OT_Rep1_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/noFAIMS_OT/OGlycoTM_EThcD_Search/noFAIMS_OT_EThcD_Rep1_bonafide.csv'
)
noFAIMS_OT_Rep2_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/noFAIMS_OT/OGlycoTM_EThcD_Search/noFAIMS_OT_EThcD_Rep2_bonafide.csv'
)
noFAIMS_OT_Rep3_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/noFAIMS_OT/OGlycoTM_EThcD_Search/noFAIMS_OT_EThcD_Rep3_bonafide.csv'
)

# noFAIMS_OT_IT
noFAIMS_OT_IT_Rep1_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/noFAIMS_OT_IT/OGlycoTM_PairedScan_Search/noFAIMS_OT_IT_Rep1_bonafide.csv'
)
noFAIMS_OT_IT_Rep2_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/noFAIMS_OT_IT/OGlycoTM_PairedScan_Search/noFAIMS_OT_IT_Rep2_bonafide.csv'
)
noFAIMS_OT_IT_Rep3_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/noFAIMS_OT_IT/OGlycoTM_PairedScan_Search/noFAIMS_OT_IT_Rep3_bonafide.csv'
)

## check glycoPSM/glycosite and glycoprotein level results
## O-GlcNAc
# FAIMS45
# glycoPSM
FAIMS45_Rep1_glycoPSM <- FAIMS45_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS45_Rep2_glycoPSM <- FAIMS45_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS45_Rep3_glycoPSM <- FAIMS45_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

# glycosite
FAIMS45_Rep1_glycosite <- FAIMS45_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS45_Rep2_glycosite <- FAIMS45_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS45_Rep3_glycosite <- FAIMS45_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

# glycoprotein
FAIMS45_Rep1_glycoprotein <- FAIMS45_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS45_Rep2_glycoprotein <- FAIMS45_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS45_Rep3_glycoprotein <- FAIMS45_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

# FAIMS50
# glycoPSM
FAIMS50_Rep1_glycoPSM <- FAIMS50_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS50_Rep2_glycoPSM <- FAIMS50_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS50_Rep3_glycoPSM <- FAIMS50_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

# glycosite
FAIMS50_Rep1_glycosite <- FAIMS50_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS50_Rep2_glycosite <- FAIMS50_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS50_Rep3_glycosite <- FAIMS50_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

# glycoprotein
FAIMS50_Rep1_glycoprotein <- FAIMS50_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS50_Rep2_glycoprotein <- FAIMS50_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS50_Rep3_glycoprotein <- FAIMS50_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

# FAIMS4045
# glycoPSM
FAIMS4045_Rep1_glycoPSM <- FAIMS4045_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS4045_Rep2_glycoPSM <- FAIMS4045_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS4045_Rep3_glycoPSM <- FAIMS4045_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

# glycosite
FAIMS4045_Rep1_glycosite <- FAIMS4045_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS4045_Rep2_glycosite <- FAIMS4045_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS4045_Rep3_glycosite <- FAIMS4045_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

# glycoprotein
FAIMS4045_Rep1_glycoprotein <- FAIMS4045_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS4045_Rep2_glycoprotein <- FAIMS4045_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS4045_Rep3_glycoprotein <- FAIMS4045_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

# FAIMS4050
# glycoPSM
FAIMS4050_Rep1_glycoPSM <- FAIMS4050_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS4050_Rep2_glycoPSM <- FAIMS4050_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS4050_Rep3_glycoPSM <- FAIMS4050_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

# glycosite
FAIMS4050_Rep1_glycosite <- FAIMS4050_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS4050_Rep2_glycosite <- FAIMS4050_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS4050_Rep3_glycosite <- FAIMS4050_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

# glycoprotein
FAIMS4050_Rep1_glycoprotein <- FAIMS4050_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS4050_Rep2_glycoprotein <- FAIMS4050_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS4050_Rep3_glycoprotein <- FAIMS4050_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

# FAIMS404550
# glycoPSM
FAIMS404550_Rep1_glycoPSM <- FAIMS404550_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS404550_Rep2_glycoPSM <- FAIMS404550_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS404550_Rep3_glycoPSM <- FAIMS404550_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

# glycosite
FAIMS404550_Rep1_glycosite <- FAIMS404550_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS404550_Rep2_glycosite <- FAIMS404550_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS404550_Rep3_glycosite <- FAIMS404550_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

# glycoprotein
FAIMS404550_Rep1_glycoprotein <- FAIMS404550_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS404550_Rep2_glycoprotein <- FAIMS404550_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS404550_Rep3_glycoprotein <- FAIMS404550_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

# FAIMS455065
# glycoPSM
FAIMS455065_Rep1_glycoPSM <- FAIMS455065_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS455065_Rep2_glycoPSM <- FAIMS455065_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS455065_Rep3_glycoPSM <- FAIMS455065_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

# glycosite
FAIMS455065_Rep1_glycosite <- FAIMS455065_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS455065_Rep2_glycosite <- FAIMS455065_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS455065_Rep3_glycosite <- FAIMS455065_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

# glycoprotein
FAIMS455065_Rep1_glycoprotein <- FAIMS455065_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS455065_Rep2_glycoprotein <- FAIMS455065_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS455065_Rep3_glycoprotein <- FAIMS455065_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

# noFAIMS_OT
# glycoPSM
noFAIMS_OT_Rep1_glycoPSM <- noFAIMS_OT_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

noFAIMS_OT_Rep2_glycoPSM <- noFAIMS_OT_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

noFAIMS_OT_Rep3_glycoPSM <- noFAIMS_OT_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

# glycosite
noFAIMS_OT_Rep1_glycosite <- noFAIMS_OT_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

noFAIMS_OT_Rep2_glycosite <- noFAIMS_OT_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

noFAIMS_OT_Rep3_glycosite <- noFAIMS_OT_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

# glycoprotein
noFAIMS_OT_Rep1_glycoprotein <- noFAIMS_OT_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

noFAIMS_OT_Rep2_glycoprotein <- noFAIMS_OT_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

noFAIMS_OT_Rep3_glycoprotein <- noFAIMS_OT_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

# noFAIMS_OT_IT
# glycoPSM
noFAIMS_OT_IT_Rep1_glycoPSM <- noFAIMS_OT_IT_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

noFAIMS_OT_IT_Rep2_glycoPSM <- noFAIMS_OT_IT_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

noFAIMS_OT_IT_Rep3_glycoPSM <- noFAIMS_OT_IT_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

# glycosite
noFAIMS_OT_IT_Rep1_glycosite <- noFAIMS_OT_IT_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

noFAIMS_OT_IT_Rep2_glycosite <- noFAIMS_OT_IT_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

noFAIMS_OT_IT_Rep3_glycosite <- noFAIMS_OT_IT_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

# glycoprotein
noFAIMS_OT_IT_Rep1_glycoprotein <- noFAIMS_OT_IT_Rep1_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

noFAIMS_OT_IT_Rep2_glycoprotein <- noFAIMS_OT_IT_Rep2_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

noFAIMS_OT_IT_Rep3_glycoprotein <- noFAIMS_OT_IT_Rep3_bonafide |> 
  filter(Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

## Mucin-Type O-Glyco
# FAIMS45
# glycoPSM
FAIMS45_Rep1_glycoPSM_Mucin <- FAIMS45_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS45_Rep2_glycoPSM_Mucin <- FAIMS45_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS45_Rep3_glycoPSM_Mucin <- FAIMS45_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

# glycosite
FAIMS45_Rep1_glycosite_Mucin <- FAIMS45_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS45_Rep2_glycosite_Mucin <- FAIMS45_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS45_Rep3_glycosite_Mucin <- FAIMS45_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

# glycoprotein
FAIMS45_Rep1_glycoprotein_Mucin <- FAIMS45_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS45_Rep2_glycoprotein_Mucin <- FAIMS45_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS45_Rep3_glycoprotein_Mucin <- FAIMS45_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

# FAIMS50
# glycoPSM
FAIMS50_Rep1_glycoPSM_Mucin <- FAIMS50_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS50_Rep2_glycoPSM_Mucin <- FAIMS50_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS50_Rep3_glycoPSM_Mucin <- FAIMS50_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

# glycosite
FAIMS50_Rep1_glycosite_Mucin <- FAIMS50_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS50_Rep2_glycosite_Mucin <- FAIMS50_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS50_Rep3_glycosite_Mucin <- FAIMS50_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

# glycoprotein
FAIMS50_Rep1_glycoprotein_Mucin <- FAIMS50_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS50_Rep2_glycoprotein_Mucin <- FAIMS50_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS50_Rep3_glycoprotein_Mucin <- FAIMS50_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

# FAIMS4045
# glycoPSM
FAIMS4045_Rep1_glycoPSM_Mucin <- FAIMS4045_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS4045_Rep2_glycoPSM_Mucin <- FAIMS4045_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS4045_Rep3_glycoPSM_Mucin <- FAIMS4045_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

# glycosite
FAIMS4045_Rep1_glycosite_Mucin <- FAIMS4045_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS4045_Rep2_glycosite_Mucin <- FAIMS4045_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS4045_Rep3_glycosite_Mucin <- FAIMS4045_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

# glycoprotein
FAIMS4045_Rep1_glycoprotein_Mucin <- FAIMS4045_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS4045_Rep2_glycoprotein_Mucin <- FAIMS4045_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS4045_Rep3_glycoprotein_Mucin <- FAIMS4045_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

# FAIMS4050
# glycoPSM
FAIMS4050_Rep1_glycoPSM_Mucin <- FAIMS4050_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS4050_Rep2_glycoPSM_Mucin <- FAIMS4050_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS4050_Rep3_glycoPSM_Mucin <- FAIMS4050_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

# glycosite
FAIMS4050_Rep1_glycosite_Mucin <- FAIMS4050_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS4050_Rep2_glycosite_Mucin <- FAIMS4050_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS4050_Rep3_glycosite_Mucin <- FAIMS4050_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

# glycoprotein
FAIMS4050_Rep1_glycoprotein_Mucin <- FAIMS4050_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS4050_Rep2_glycoprotein_Mucin <- FAIMS4050_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS4050_Rep3_glycoprotein_Mucin <- FAIMS4050_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

# FAIMS404550
# glycoPSM
FAIMS404550_Rep1_glycoPSM_Mucin <- FAIMS404550_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS404550_Rep2_glycoPSM_Mucin <- FAIMS404550_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS404550_Rep3_glycoPSM_Mucin <- FAIMS404550_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

# glycosite
FAIMS404550_Rep1_glycosite_Mucin <- FAIMS404550_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS404550_Rep2_glycosite_Mucin <- FAIMS404550_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS404550_Rep3_glycosite_Mucin <- FAIMS404550_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

# glycoprotein
FAIMS404550_Rep1_glycoprotein_Mucin <- FAIMS404550_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS404550_Rep2_glycoprotein_Mucin <- FAIMS404550_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS404550_Rep3_glycoprotein_Mucin <- FAIMS404550_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

# FAIMS455065
# glycoPSM
FAIMS455065_Rep1_glycoPSM_Mucin <- FAIMS455065_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS455065_Rep2_glycoPSM_Mucin <- FAIMS455065_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

FAIMS455065_Rep3_glycoPSM_Mucin <- FAIMS455065_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

# glycosite
FAIMS455065_Rep1_glycosite_Mucin <- FAIMS455065_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS455065_Rep2_glycosite_Mucin <- FAIMS455065_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

FAIMS455065_Rep3_glycosite_Mucin <- FAIMS455065_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

# glycoprotein
FAIMS455065_Rep1_glycoprotein_Mucin <- FAIMS455065_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS455065_Rep2_glycoprotein_Mucin <- FAIMS455065_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

FAIMS455065_Rep3_glycoprotein_Mucin <- FAIMS455065_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

# noFAIMS_OT
# glycoPSM
noFAIMS_OT_Rep1_glycoPSM_Mucin <- noFAIMS_OT_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

noFAIMS_OT_Rep2_glycoPSM_Mucin <- noFAIMS_OT_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

noFAIMS_OT_Rep3_glycoPSM_Mucin <- noFAIMS_OT_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

# glycosite
noFAIMS_OT_Rep1_glycosite_Mucin <- noFAIMS_OT_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

noFAIMS_OT_Rep2_glycosite_Mucin <- noFAIMS_OT_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

noFAIMS_OT_Rep3_glycosite_Mucin <- noFAIMS_OT_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

# glycoprotein
noFAIMS_OT_Rep1_glycoprotein_Mucin <- noFAIMS_OT_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

noFAIMS_OT_Rep2_glycoprotein_Mucin <- noFAIMS_OT_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

noFAIMS_OT_Rep3_glycoprotein_Mucin <- noFAIMS_OT_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

# noFAIMS_OT_IT
# glycoPSM
noFAIMS_OT_IT_Rep1_glycoPSM_Mucin <- noFAIMS_OT_IT_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

noFAIMS_OT_IT_Rep2_glycoPSM_Mucin <- noFAIMS_OT_IT_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

noFAIMS_OT_IT_Rep3_glycoPSM_Mucin <- noFAIMS_OT_IT_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count()

# glycosite
noFAIMS_OT_IT_Rep1_glycosite_Mucin <- noFAIMS_OT_IT_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

noFAIMS_OT_IT_Rep2_glycosite_Mucin <- noFAIMS_OT_IT_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

noFAIMS_OT_IT_Rep3_glycosite_Mucin <- noFAIMS_OT_IT_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  count(Confidence.Level)

# glycoprotein
noFAIMS_OT_IT_Rep1_glycoprotein_Mucin <- noFAIMS_OT_IT_Rep1_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

noFAIMS_OT_IT_Rep2_glycoprotein_Mucin <- noFAIMS_OT_IT_Rep2_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

noFAIMS_OT_IT_Rep3_glycoprotein_Mucin <- noFAIMS_OT_IT_Rep3_bonafide |> 
  filter(!Total.Glycan.Composition == 'HexNAt(1)') |> 
  distinct(Protein.ID) |> 
  nrow()

## barplot data
# O-GlcNAc data
# glycoPSM data
glycoPSM_OGlcNAc <- tibble(
  Condition = rep(c("FAIMS45", "FAIMS50", "FAIMS4045", "FAIMS4050", 
                    "FAIMS404550", "FAIMS455065", "noFAIMS_OT", "noFAIMS_OT_IT"), each = 3),
  Replicate = rep(c("Rep1", "Rep2", "Rep3"), 8),
  Count = c(
    FAIMS45_Rep1_glycoPSM$n, FAIMS45_Rep2_glycoPSM$n, FAIMS45_Rep3_glycoPSM$n,
    FAIMS50_Rep1_glycoPSM$n, FAIMS50_Rep2_glycoPSM$n, FAIMS50_Rep3_glycoPSM$n,
    FAIMS4045_Rep1_glycoPSM$n, FAIMS4045_Rep2_glycoPSM$n, FAIMS4045_Rep3_glycoPSM$n,
    FAIMS4050_Rep1_glycoPSM$n, FAIMS4050_Rep2_glycoPSM$n, FAIMS4050_Rep3_glycoPSM$n,
    FAIMS404550_Rep1_glycoPSM$n, FAIMS404550_Rep2_glycoPSM$n, FAIMS404550_Rep3_glycoPSM$n,
    FAIMS455065_Rep1_glycoPSM$n, FAIMS455065_Rep2_glycoPSM$n, FAIMS455065_Rep3_glycoPSM$n,
    noFAIMS_OT_Rep1_glycoPSM$n, noFAIMS_OT_Rep2_glycoPSM$n, noFAIMS_OT_Rep3_glycoPSM$n,
    noFAIMS_OT_IT_Rep1_glycoPSM$n, noFAIMS_OT_IT_Rep2_glycoPSM$n, noFAIMS_OT_IT_Rep3_glycoPSM$n
  )
)

glycoPSM_OGlcNAc_mean <- glycoPSM_OGlcNAc %>%
  group_by(Condition) %>%
  summarise(Mean = mean(Count))

# glycosite data
glycosite_OGlcNAc <- bind_rows(
  FAIMS45_Rep1_glycosite %>% mutate(Condition = "FAIMS45", Replicate = "Rep1"),
  FAIMS45_Rep2_glycosite %>% mutate(Condition = "FAIMS45", Replicate = "Rep2"),
  FAIMS45_Rep3_glycosite %>% mutate(Condition = "FAIMS45", Replicate = "Rep3"),
  FAIMS50_Rep1_glycosite %>% mutate(Condition = "FAIMS50", Replicate = "Rep1"),
  FAIMS50_Rep2_glycosite %>% mutate(Condition = "FAIMS50", Replicate = "Rep2"),
  FAIMS50_Rep3_glycosite %>% mutate(Condition = "FAIMS50", Replicate = "Rep3"),
  FAIMS4045_Rep1_glycosite %>% mutate(Condition = "FAIMS4045", Replicate = "Rep1"),
  FAIMS4045_Rep2_glycosite %>% mutate(Condition = "FAIMS4045", Replicate = "Rep2"),
  FAIMS4045_Rep3_glycosite %>% mutate(Condition = "FAIMS4045", Replicate = "Rep3"),
  FAIMS4050_Rep1_glycosite %>% mutate(Condition = "FAIMS4050", Replicate = "Rep1"),
  FAIMS4050_Rep2_glycosite %>% mutate(Condition = "FAIMS4050", Replicate = "Rep2"),
  FAIMS4050_Rep3_glycosite %>% mutate(Condition = "FAIMS4050", Replicate = "Rep3"),
  FAIMS404550_Rep1_glycosite %>% mutate(Condition = "FAIMS404550", Replicate = "Rep1"),
  FAIMS404550_Rep2_glycosite %>% mutate(Condition = "FAIMS404550", Replicate = "Rep2"),
  FAIMS404550_Rep3_glycosite %>% mutate(Condition = "FAIMS404550", Replicate = "Rep3"),
  FAIMS455065_Rep1_glycosite %>% mutate(Condition = "FAIMS455065", Replicate = "Rep1"),
  FAIMS455065_Rep2_glycosite %>% mutate(Condition = "FAIMS455065", Replicate = "Rep2"),
  FAIMS455065_Rep3_glycosite %>% mutate(Condition = "FAIMS455065", Replicate = "Rep3"),
  noFAIMS_OT_Rep1_glycosite %>% mutate(Condition = "noFAIMS_OT", Replicate = "Rep1"),
  noFAIMS_OT_Rep2_glycosite %>% mutate(Condition = "noFAIMS_OT", Replicate = "Rep2"),
  noFAIMS_OT_Rep3_glycosite %>% mutate(Condition = "noFAIMS_OT", Replicate = "Rep3"),
  noFAIMS_OT_IT_Rep1_glycosite %>% mutate(Condition = "noFAIMS_OT_IT", Replicate = "Rep1"),
  noFAIMS_OT_IT_Rep2_glycosite %>% mutate(Condition = "noFAIMS_OT_IT", Replicate = "Rep2"),
  noFAIMS_OT_IT_Rep3_glycosite %>% mutate(Condition = "noFAIMS_OT_IT", Replicate = "Rep3")
) %>%
  rename(Count = n, Confidence_Level = Confidence.Level)

glycosite_OGlcNAc_mean <- glycosite_OGlcNAc %>%
  group_by(Condition, Confidence_Level) %>%
  summarise(Mean = mean(Count), .groups = "drop")

# glycoprotein data
glycoprotein_OGlcNAc <- tibble(
  Condition = rep(c("FAIMS45", "FAIMS50", "FAIMS4045", "FAIMS4050", 
                    "FAIMS404550", "FAIMS455065", "noFAIMS_OT", "noFAIMS_OT_IT"), each = 3),
  Replicate = rep(c("Rep1", "Rep2", "Rep3"), 8),
  Count = c(
    FAIMS45_Rep1_glycoprotein, FAIMS45_Rep2_glycoprotein, FAIMS45_Rep3_glycoprotein,
    FAIMS50_Rep1_glycoprotein, FAIMS50_Rep2_glycoprotein, FAIMS50_Rep3_glycoprotein,
    FAIMS4045_Rep1_glycoprotein, FAIMS4045_Rep2_glycoprotein, FAIMS4045_Rep3_glycoprotein,
    FAIMS4050_Rep1_glycoprotein, FAIMS4050_Rep2_glycoprotein, FAIMS4050_Rep3_glycoprotein,
    FAIMS404550_Rep1_glycoprotein, FAIMS404550_Rep2_glycoprotein, FAIMS404550_Rep3_glycoprotein,
    FAIMS455065_Rep1_glycoprotein, FAIMS455065_Rep2_glycoprotein, FAIMS455065_Rep3_glycoprotein,
    noFAIMS_OT_Rep1_glycoprotein, noFAIMS_OT_Rep2_glycoprotein, noFAIMS_OT_Rep3_glycoprotein,
    noFAIMS_OT_IT_Rep1_glycoprotein, noFAIMS_OT_IT_Rep2_glycoprotein, noFAIMS_OT_IT_Rep3_glycoprotein
  )
)

glycoprotein_OGlcNAc_mean <- glycoprotein_OGlcNAc %>%
  group_by(Condition) %>%
  summarise(Mean = mean(Count))

# Mucin-Type O-Glyco data
# glycoPSM data
glycoPSM_Mucin <- tibble(
  Condition = rep(c("FAIMS45", "FAIMS50", "FAIMS4045", "FAIMS4050", 
                    "FAIMS404550", "FAIMS455065", "noFAIMS_OT", "noFAIMS_OT_IT"), each = 3),
  Replicate = rep(c("Rep1", "Rep2", "Rep3"), 8),
  Count = c(
    FAIMS45_Rep1_glycoPSM_Mucin$n, FAIMS45_Rep2_glycoPSM_Mucin$n, FAIMS45_Rep3_glycoPSM_Mucin$n,
    FAIMS50_Rep1_glycoPSM_Mucin$n, FAIMS50_Rep2_glycoPSM_Mucin$n, FAIMS50_Rep3_glycoPSM_Mucin$n,
    FAIMS4045_Rep1_glycoPSM_Mucin$n, FAIMS4045_Rep2_glycoPSM_Mucin$n, FAIMS4045_Rep3_glycoPSM_Mucin$n,
    FAIMS4050_Rep1_glycoPSM_Mucin$n, FAIMS4050_Rep2_glycoPSM_Mucin$n, FAIMS4050_Rep3_glycoPSM_Mucin$n,
    FAIMS404550_Rep1_glycoPSM_Mucin$n, FAIMS404550_Rep2_glycoPSM_Mucin$n, FAIMS404550_Rep3_glycoPSM_Mucin$n,
    FAIMS455065_Rep1_glycoPSM_Mucin$n, FAIMS455065_Rep2_glycoPSM_Mucin$n, FAIMS455065_Rep3_glycoPSM_Mucin$n,
    noFAIMS_OT_Rep1_glycoPSM_Mucin$n, noFAIMS_OT_Rep2_glycoPSM_Mucin$n, noFAIMS_OT_Rep3_glycoPSM_Mucin$n,
    noFAIMS_OT_IT_Rep1_glycoPSM_Mucin$n, noFAIMS_OT_IT_Rep2_glycoPSM_Mucin$n, noFAIMS_OT_IT_Rep3_glycoPSM_Mucin$n
  )
)

glycoPSM_Mucin_mean <- glycoPSM_Mucin %>%
  group_by(Condition) %>%
  summarise(Mean = mean(Count))

# glycosite data
glycosite_Mucin <- bind_rows(
  FAIMS45_Rep1_glycosite_Mucin %>% mutate(Condition = "FAIMS45", Replicate = "Rep1"),
  FAIMS45_Rep2_glycosite_Mucin %>% mutate(Condition = "FAIMS45", Replicate = "Rep2"),
  FAIMS45_Rep3_glycosite_Mucin %>% mutate(Condition = "FAIMS45", Replicate = "Rep3"),
  FAIMS50_Rep1_glycosite_Mucin %>% mutate(Condition = "FAIMS50", Replicate = "Rep1"),
  FAIMS50_Rep2_glycosite_Mucin %>% mutate(Condition = "FAIMS50", Replicate = "Rep2"),
  FAIMS50_Rep3_glycosite_Mucin %>% mutate(Condition = "FAIMS50", Replicate = "Rep3"),
  FAIMS4045_Rep1_glycosite_Mucin %>% mutate(Condition = "FAIMS4045", Replicate = "Rep1"),
  FAIMS4045_Rep2_glycosite_Mucin %>% mutate(Condition = "FAIMS4045", Replicate = "Rep2"),
  FAIMS4045_Rep3_glycosite_Mucin %>% mutate(Condition = "FAIMS4045", Replicate = "Rep3"),
  FAIMS4050_Rep1_glycosite_Mucin %>% mutate(Condition = "FAIMS4050", Replicate = "Rep1"),
  FAIMS4050_Rep2_glycosite_Mucin %>% mutate(Condition = "FAIMS4050", Replicate = "Rep2"),
  FAIMS4050_Rep3_glycosite_Mucin %>% mutate(Condition = "FAIMS4050", Replicate = "Rep3"),
  FAIMS404550_Rep1_glycosite_Mucin %>% mutate(Condition = "FAIMS404550", Replicate = "Rep1"),
  FAIMS404550_Rep2_glycosite_Mucin %>% mutate(Condition = "FAIMS404550", Replicate = "Rep2"),
  FAIMS404550_Rep3_glycosite_Mucin %>% mutate(Condition = "FAIMS404550", Replicate = "Rep3"),
  FAIMS455065_Rep1_glycosite_Mucin %>% mutate(Condition = "FAIMS455065", Replicate = "Rep1"),
  FAIMS455065_Rep2_glycosite_Mucin %>% mutate(Condition = "FAIMS455065", Replicate = "Rep2"),
  FAIMS455065_Rep3_glycosite_Mucin %>% mutate(Condition = "FAIMS455065", Replicate = "Rep3"),
  noFAIMS_OT_Rep1_glycosite_Mucin %>% mutate(Condition = "noFAIMS_OT", Replicate = "Rep1"),
  noFAIMS_OT_Rep2_glycosite_Mucin %>% mutate(Condition = "noFAIMS_OT", Replicate = "Rep2"),
  noFAIMS_OT_Rep3_glycosite_Mucin %>% mutate(Condition = "noFAIMS_OT", Replicate = "Rep3"),
  noFAIMS_OT_IT_Rep1_glycosite_Mucin %>% mutate(Condition = "noFAIMS_OT_IT", Replicate = "Rep1"),
  noFAIMS_OT_IT_Rep2_glycosite_Mucin %>% mutate(Condition = "noFAIMS_OT_IT", Replicate = "Rep2"),
  noFAIMS_OT_IT_Rep3_glycosite_Mucin %>% mutate(Condition = "noFAIMS_OT_IT", Replicate = "Rep3")
) %>%
  rename(Count = n, Confidence_Level = Confidence.Level)

glycosite_Mucin_mean <- glycosite_Mucin %>%
  group_by(Condition, Confidence_Level) %>%
  summarise(Mean = mean(Count), .groups = "drop")

# glycoprotein data
glycoprotein_Mucin <- tibble(
  Condition = rep(c("FAIMS45", "FAIMS50", "FAIMS4045", "FAIMS4050", 
                    "FAIMS404550", "FAIMS455065", "noFAIMS_OT", "noFAIMS_OT_IT"), each = 3),
  Replicate = rep(c("Rep1", "Rep2", "Rep3"), 8),
  Count = c(
    FAIMS45_Rep1_glycoprotein_Mucin, FAIMS45_Rep2_glycoprotein_Mucin, FAIMS45_Rep3_glycoprotein_Mucin,
    FAIMS50_Rep1_glycoprotein_Mucin, FAIMS50_Rep2_glycoprotein_Mucin, FAIMS50_Rep3_glycoprotein_Mucin,
    FAIMS4045_Rep1_glycoprotein_Mucin, FAIMS4045_Rep2_glycoprotein_Mucin, FAIMS4045_Rep3_glycoprotein_Mucin,
    FAIMS4050_Rep1_glycoprotein_Mucin, FAIMS4050_Rep2_glycoprotein_Mucin, FAIMS4050_Rep3_glycoprotein_Mucin,
    FAIMS404550_Rep1_glycoprotein_Mucin, FAIMS404550_Rep2_glycoprotein_Mucin, FAIMS404550_Rep3_glycoprotein_Mucin,
    FAIMS455065_Rep1_glycoprotein_Mucin, FAIMS455065_Rep2_glycoprotein_Mucin, FAIMS455065_Rep3_glycoprotein_Mucin,
    noFAIMS_OT_Rep1_glycoprotein_Mucin, noFAIMS_OT_Rep2_glycoprotein_Mucin, noFAIMS_OT_Rep3_glycoprotein_Mucin,
    noFAIMS_OT_IT_Rep1_glycoprotein_Mucin, noFAIMS_OT_IT_Rep2_glycoprotein_Mucin, noFAIMS_OT_IT_Rep3_glycoprotein_Mucin
  )
)

glycoprotein_Mucin_mean <- glycoprotein_Mucin %>%
  group_by(Condition) %>%
  summarise(Mean = mean(Count))

## barplot
# O-GlcNAc glycoPSM
OGlcNAc_glycoPSM <- ggplot(glycoPSM_OGlcNAc_mean, aes(x = Condition, y = Mean)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_point(data = glycoPSM_OGlcNAc, aes(x = Condition, y = Count), 
             size = 1, color = "black") +
  labs(title = "O-GlcNAc: Total glycoPSM", x = "", y = "Count") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1, color = 'black'),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 8, color = 'black'),
        title = element_text(size = 8, color = 'black'))

ggsave(
  filename = '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/OGlcNAc_glycoPSM.tif',
  plot = OGlcNAc_glycoPSM,
  height = 3, width = 3, 
  dpi = 300
)

# O-GlcNAc glycosite
OGlcNAc_glycosite <- ggplot(glycosite_OGlcNAc_mean, aes(x = Condition, y = Mean, fill = Confidence_Level)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.7) +
  geom_point(data = glycosite_OGlcNAc, aes(x = Condition, y = Count, group = Confidence_Level),
             position = position_dodge(width = 0.9), size = 1, color = "black") +
  labs(title = "O-GlcNAc: Total glycosite", x = "", y = "Count", fill = "Confidence Level") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1, color = 'black'),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 8, color = 'black'),
        title = element_text(size = 8, color = 'black'),
        legend.text = element_text(size = 7, color = "black"),
        legend.key.size = unit(0.1, 'in'))

ggsave(
  filename = '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/OGlcNAc_glycosite.tif',
  plot = OGlcNAc_glycosite,
  height = 3, width = 4, 
  dpi = 300
)

# O-GlcNAc glycoprotein
OGlcNAc_glycoprotein <- ggplot(glycoprotein_OGlcNAc_mean, aes(x = Condition, y = Mean)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_point(data = glycoprotein_OGlcNAc, aes(x = Condition, y = Count), 
             size = 1, color = "black") +
  labs(title = "O-GlcNAc: Unique glycoprotein", x = "", y = "Count") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1, color = 'black'),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 8, color = 'black'),
        title = element_text(size = 8, color = 'black'))

ggsave(
  filename = '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/OGlcNAc_glycoprotein.tif',
  plot = OGlcNAc_glycoprotein,
  height = 3, width = 3, 
  dpi = 300
)

# Mucin-Type O-Glyco glycoPSM
Mucin_glycoPSM <- ggplot(glycoPSM_Mucin_mean, aes(x = Condition, y = Mean)) +
  geom_bar(stat = "identity", fill = "coral", alpha = 0.7) +
  geom_point(data = glycoPSM_Mucin, aes(x = Condition, y = Count), 
             size = 1, color = "black") +
  labs(title = "Mucin-Type O-Glyco: Total glycoPSM", x = "", y = "Count") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1, color = 'black'),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 8, color = 'black'),
        title = element_text(size = 8, color = 'black'))

ggsave(
  filename = '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/Mucin_glycoPSM.tif',
  plot = Mucin_glycoPSM,
  height = 3, width = 3, 
  dpi = 300
)

# Mucin-Type O-Glyco glycosite
Mucin_glycosite <- ggplot(glycosite_Mucin_mean, aes(x = Condition, y = Mean, fill = Confidence_Level)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.7) +
  geom_point(data = glycosite_Mucin, aes(x = Condition, y = Count, group = Confidence_Level),
             position = position_dodge(width = 0.9), size = 1, color = "black") +
  labs(title = "Mucin-Type O-Glyco: Total glycosite", x = "", y = "Count", fill = "Confidence Level") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1, color = 'black'),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 8, color = 'black'),
        title = element_text(size = 8, color = 'black'),
        legend.text = element_text(size = 7, color = "black"),
        legend.key.size = unit(0.1, 'in'))

ggsave(
  filename = '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/Mucin_glycosite.tif',
  plot = Mucin_glycosite,
  height = 3, width = 4, 
  dpi = 300
)

# Mucin-Type glycoprotein
Mucin_glycoprotein <- ggplot(glycoprotein_Mucin_mean, aes(x = Condition, y = Mean)) +
  geom_bar(stat = "identity", fill = "coral", alpha = 0.7) +
  geom_point(data = glycoprotein_Mucin, aes(x = Condition, y = Count), 
             size = 1, color = "black") +
  labs(title = "Mucin-Type O-Glyco: Unique glycoprotein", x = "", y = "Count") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1, color = 'black'),
        axis.text.y = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 8, color = 'black'),
        title = element_text(size = 8, color = 'black'))

ggsave(
  filename = '/Volumes/cos-lab-rwu60/Longping/HeLaGlyco_09272025/Mucin_glycoprotein.tif',
  plot = Mucin_glycoprotein,
  height = 3, width = 3, 
  dpi = 300
)
