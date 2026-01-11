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
  paste0(source_file_path, 'raw/OGlyco_HEK293T_raw.csv')
)
OGlyco_HepG2_raw <- read_csv(
  paste0(source_file_path, 'raw/OGlyco_HepG2_raw.csv')
)
OGlyco_Jurkat_raw <- read_csv(
  paste0(source_file_path, 'raw/OGlyco_Jurkat_raw.csv')
)

# bonafide glyco data
OGlyco_HEK293T_bonafide <- read_csv(
  paste0(source_file_path, 'filtered/OGlyco_HEK293T_bonafide.csv')
)
OGlyco_HepG2_bonafide <- read_csv(
  paste0(source_file_path, 'filtered/OGlyco_HepG2_bonafide.csv')
)
OGlyco_Jurkat_bonafide <- read_csv(
  paste0(source_file_path, 'filtered/OGlyco_Jurkat_bonafide.csv')
)

# O-GlcNAc filtered data
OGlcNAc_HEK293T <- read_csv(
  paste0(source_file_path, 'filtered/OGlcNAc_HEK293T.csv')
)
OGlcNAc_HepG2 <- read_csv(
  paste0(source_file_path, 'filtered/OGlcNAc_HepG2.csv')
)
OGlcNAc_Jurkat <- read_csv(
  paste0(source_file_path, 'filtered/OGlcNAc_Jurkat.csv')
)

# localized OGlyco site
OGlyco_site_HEK293T <- read_csv(
  paste0(source_file_path, 'site/OGlyco_site_HEK293T.csv')
)
OGlyco_site_HepG2 <- read_csv(
  paste0(source_file_path, 'site/OGlyco_site_HepG2.csv')
)
OGlyco_site_Jurkat <- read_csv(
  paste0(source_file_path, 'site/OGlyco_site_Jurkat.csv')
)

# O-GlcNAc site filtered data
OGlcNAc_site_HEK293T <- read_csv(
  paste0(source_file_path, 'site/OGlcNAc_site_HEK293T.csv')
)
OGlcNAc_site_HepG2 <- read_csv(
  paste0(source_file_path, 'site/OGlcNAc_site_HepG2.csv')
)
OGlcNAc_site_Jurkat <- read_csv(
  paste0(source_file_path, 'site/OGlcNAc_site_Jurkat.csv')
)

# O-GlcNAc protein level quantification data
OGlcNAc_protein_quant_HEK293T <- read_csv(
  paste0(source_file_path, 'quantification/OGlcNAc_protein_quant_HEK293T.csv')
)
OGlcNAc_protein_quant_HepG2 <- read_csv(
  paste0(source_file_path, 'quantification/OGlcNAc_protein_quant_HepG2.csv')
)
OGlcNAc_protein_quant_Jurkat <- read_csv(
  paste0(source_file_path, 'quantification/OGlcNAc_protein_quant_Jurkat.csv')
)

# O-GlcNAc site level quantification data
OGlcNAc_site_quant_HEK293T <- read_csv(
  paste0(source_file_path, 'quantification/OGlcNAc_site_quant_HEK293T.csv')
)
OGlcNAc_site_quant_HepG2 <- read_csv(
  paste0(source_file_path, 'quantification/OGlcNAc_site_quant_HepG2.csv')
)
OGlcNAc_site_quant_Jurkat <- read_csv(
  paste0(source_file_path, 'quantification/OGlcNAc_site_quant_Jurkat.csv')
)

# OGlcNAc protein total
OGlcNAc_protein_total <- read_csv(
  paste0(source_file_path, 'protein_lists/OGlcNAc_protein_total.csv')
)

# OGlcNAc site total
OGlcNAc_site_total <- read_csv(
  paste0(source_file_path, 'protein_lists/OGlcNAc_site_total.csv')
)

# Whole proteome protein level quantification data
WP_protein_quant_HEK293T <- read_csv(
  paste0(source_file_path, 'quantification/WP_protein_quant_HEK293T.csv')
)
WP_protein_quant_HepG2 <- read_csv(
  paste0(source_file_path, 'quantification/WP_protein_quant_HepG2.csv')
)
WP_protein_quant_Jurkat <- read_csv(
  paste0(source_file_path, 'quantification/WP_protein_quant_Jurkat.csv')
)

