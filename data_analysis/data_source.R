# import packages
library(tidyverse)

# Data roots.
#
# The network share is canonical (see 00_FILE_MAP.md at its root). It is also SMB and drops reads
# intermittently, so a curated copy of data_source/ and Figures/ is kept on the Expansion drive.
# Point OGLYCOTM_DATA / OGLYCOTM_FIGURES at that copy to read from it instead:
#
#   Sys.setenv(OGLYCOTM_DATA    = "/Volumes/Expansion/Longping/OGlycoTM/data_source")
#   Sys.setenv(OGLYCOTM_FIGURES = "/Volumes/Expansion/Longping/OGlycoTM/Figures")
#
# Same idiom as OGLYCO_DATA in export_web_data.py. Unset (or unreadable) falls back to the
# network path, so every existing script keeps working unchanged.
resolve_data_root <- function(env_var, default) {
  p <- Sys.getenv(env_var, unset = "")
  if (!nzchar(p)) return(default)
  if (!dir.exists(p)) {
    warning(env_var, " is set to '", p, "' but that directory is not readable; ",
            "falling back to ", default, call. = FALSE)
    return(default)
  }
  if (!grepl("/$", p)) p <- paste0(p, "/")
  message("data_source.R: ", env_var, " -> ", p)
  p
}

# source file path
source_file_path <- resolve_data_root(
  "OGLYCOTM_DATA",
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'
)

# figure file path
figure_file_path <- resolve_data_root(
  "OGLYCOTM_FIGURES",
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/'
)

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

