# Data Source for Differential Analysis Results
# Source this file for analyses that need DE data

# Source base data paths and colors
source('data_source.R')

# =============================================================================
# O-GlcNAc Protein Level Differential Analysis Results
# =============================================================================

OGlcNAc_protein_DE_HEK293T <- read_csv(
  paste0(source_file_path, 'differential_analysis/OGlcNAc_protein_DE_HEK293T.csv')
)
OGlcNAc_protein_DE_HepG2 <- read_csv(
  paste0(source_file_path, 'differential_analysis/OGlcNAc_protein_DE_HepG2.csv')
)
OGlcNAc_protein_DE_Jurkat <- read_csv(
  paste0(source_file_path, 'differential_analysis/OGlcNAc_protein_DE_Jurkat.csv')
)

# =============================================================================
# O-GlcNAc Site Level Differential Analysis Results
# =============================================================================

OGlcNAc_site_DE_HEK293T <- read_csv(
  paste0(source_file_path, 'differential_analysis/OGlcNAc_site_DE_HEK293T.csv')
)
OGlcNAc_site_DE_HepG2 <- read_csv(
  paste0(source_file_path, 'differential_analysis/OGlcNAc_site_DE_HepG2.csv')
)
OGlcNAc_site_DE_Jurkat <- read_csv(
  paste0(source_file_path, 'differential_analysis/OGlcNAc_site_DE_Jurkat.csv')
)

# =============================================================================
# Commonly Regulated O-GlcNAc Proteins
# =============================================================================

OGlcNAc_protein_commonly_up <- read_csv(
  paste0(source_file_path, 'differential_analysis/OGlcNAc_protein_commonly_up.csv')
)
OGlcNAc_protein_commonly_down <- read_csv(
  paste0(source_file_path, 'differential_analysis/OGlcNAc_protein_commonly_down.csv')
)

# =============================================================================
# Whole Proteome Protein Level Differential Analysis Results
# =============================================================================

WP_protein_DE_HEK293T <- read_csv(
  paste0(source_file_path, 'differential_analysis/WP_protein_DE_HEK293T.csv')
)
WP_protein_DE_HepG2 <- read_csv(
  paste0(source_file_path, 'differential_analysis/WP_protein_DE_HepG2.csv')
)
WP_protein_DE_Jurkat <- read_csv(
  paste0(source_file_path, 'differential_analysis/WP_protein_DE_Jurkat.csv')
)
