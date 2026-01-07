# import packages
library(tidyverse)

# source file path
source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'

# figure file path
figure_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/'

# define color palette
color_palette <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4")

# named version for specific uses
colors_glycan <- c("O-GlcNAc" = "#F39B7F", "O-GalNAc" = "#4DBBD5")
colors_cell <- c("HEK293T" = "#4DBBD5", "HepG2" = "#F39B7F", "Jurkat" = "#00A087")

# raw data
OGlyco_HEK293T_raw <- read_csv(
  paste0(source_file_path, 'OGlyco_HEK293T_raw.csv')
)
OGlyco_HepG2_raw <- read_csv(
  paste0(source_file_path, 'OGlyco_HepG2_raw.csv')
)
OGlyco_Jurkat_raw <- read_csv(
  paste0(source_file_path, 'OGlyco_Jurkat_raw.csv')
)

# bonafide glyco data
OGlyco_HEK293T_bonafide <- read_csv(
  paste0(source_file_path, 'OGlyco_HEK293T_bonafide.csv')
)
OGlyco_HepG2_bonafide <- read_csv(
  paste0(source_file_path, 'OGlyco_HepG2_bonafide.csv')
)
OGlyco_Jurkat_bonafide <- read_csv(
  paste0(source_file_path, 'OGlyco_Jurkat_bonafide.csv')
)

# localized OGlyco site
OGlyco_site_HEK293T <- read_csv(
  paste0(source_file_path, 'OGlyco_site_HEK293T.csv')
)
OGlyco_site_HepG2 <- read_csv(
  paste0(source_file_path, 'OGlyco_site_HepG2.csv')
)
OGlyco_site_Jurkat <- read_csv(
  paste0(source_file_path, 'OGlyco_site_Jurkat.csv')
)
