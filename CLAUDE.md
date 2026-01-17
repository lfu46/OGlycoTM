# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

OGlycoTM is an R-based data analysis project for analyzing O-glycosylation (O-GlcNAc and O-GalNAc modifications) in mass spectrometry proteomics data across three human cell types: HEK293T, HepG2, and Jurkat. This is a research manuscript project generating publication-quality figures and statistical analysis.

## Data Pipeline Architecture

The analysis follows a sequential pipeline with two parallel tracks:

### O-Glycoproteomics Pipeline
1. **data_import.R** - Imports raw TSV files from OGlycoTM mass spectrometry search results, filters for human proteins, removes decoys/contaminants
2. **data_filtering.R** - Filters for high-confidence localized sites (Level1, Level1b), removes cysteine artifacts, creates "bonafide" datasets
3. **data_quantification.R** - Aggregates PSM intensities to protein and site levels for O-GlcNAc
4. **data_normalization.R** - Applies SL (sample loading) + TMM normalization using edgeR
5. **differential_analysis.R** - Performs limma-based differential expression (Tuni vs Ctrl)

### Whole Proteome Pipeline
1. **data_import.R** - Imports WP raw CSV files, renames TMT channels (Sn.126-131)
2. **data_filtering.R** - Quality filters: XCorr > 1.2, PPM -10 to 10, all Sn > 5, removes decoys/contaminants
3. **data_quantification.R** - Aggregates to protein level by UniProt_Accession
4. **data_normalization.R** - Same SL + TMM normalization
5. **differential_analysis.R** - Same limma analysis

### Data Source Files
- **data_source.R** - Central configuration with file paths, color palettes, loads all filtered/quantified datasets
- **data_source_DE.R** - Loads differential analysis results (sources data_source.R first)

## Running Scripts

```r
# Set working directory to data_analysis folder
setwd("data_analysis")

# For analysis using pre-processed data:
source('data_source.R')           # Load base data and config
source('data_source_DE.R')        # Load differential analysis results
source('Figure2.R')               # Generate Figure 2 panels
source('Figure3.R')               # Generate Figure 3 panels

# For full pipeline from raw data (run in order):
source('data_import.R')           # Import raw data
source('data_filtering.R')        # Filter & curate
source('data_quantification.R')   # Aggregate to protein/site level
source('data_normalization.R')    # Normalize intensities
source('differential_analysis.R') # Differential expression analysis
```

## Key Configuration (data_source.R)

- **source_file_path**: `/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/`
- **figure_file_path**: `/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/`
- **Color palettes**:
  - `colors_glycan`: O-GlcNAc (#F39B7F salmon), O-GalNAc (#4DBBD5 blue)
  - `colors_cell`: HEK293T (#4DBBD5), HepG2 (#F39B7F), Jurkat (#00A087)

## Experimental Design

- **Conditions**: Tuni (tunicamycin treatment) vs Ctrl (control)
- **Replicates**: 3 per condition (channels 126-128 = Tuni, 129-131 = Ctrl)
- **TMT channels**: Intensity.Tuni_1, Tuni_2, Tuni_3, Ctrl_4, Ctrl_5, Ctrl_6
- **Normalized columns**: `_sl` suffix for SL-normalized, `_sl_tmm` suffix for SL+TMM normalized

## Key R Packages

- `tidyverse` - Data manipulation and ggplot2 visualization
- `edgeR` - TMM normalization (calcNormFactors)
- `limma` - Differential expression analysis
- `clusterProfiler` - GO enrichment analysis
- `eulerr` - Proportional Venn/Euler diagrams
- `introdataviz` - Split violin plots (install from GitHub: `remotes::install_github("psyteachr/introdataviz")`)
- `ggpubr`, `rstatix` - Statistical annotations on plots

## Data Files

All data files (CSV, TSV, XLSX) are in .gitignore. Data is stored on external network drive at `/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/`.

**Key data subfolders:**
- `raw/` - Raw imported data
- `filtered/` - Quality-filtered datasets
- `quantification/` - Protein/site level quantified data
- `normalization/` - Normalized intensity data
- `differential_analysis/` - DE results and commonly regulated protein lists
- `enrichment/` - GO enrichment results

## Figure Scripts

Each figure script can run independently by sourcing `data_source.R` or `data_source_DE.R` first:
- **Figure2.R** - O-GlcNAc identification overview (Euler diagrams, site distribution)
- **Figure3.R** - Differential expression analysis (volcano plots, heatmaps)
- **Figure4.R** - Cell-type specific analysis (circular heatmaps, GO enrichment)
- **Figure5.R** - Subcellular localization analysis (proportion barplots, location-specific dotplots)
- **Figure6.R** - Structural feature analysis (logFC distribution, secondary structure, IDR effects, site-specific examples)

## PyMOL 3D Structure Visualization (Figure 6F)

PyMOL scripts for protein structure visualization with O-GlcNAc site highlighting:
- **Figure6F_pymol.py** - Main Python script for PyMOL, colors by pLDDT, highlights sites
- **Figure6F_ray_modes.py** - Comparison of ray trace modes (0-3)

```bash
# Run PyMOL script (requires PyMOL installed via Homebrew)
/opt/homebrew/bin/pymol -c -q Figure6F_pymol.py

# AlphaFold structures are downloaded from:
# https://alphafold.ebi.ac.uk/files/AF-{UniProt_ID}-F1-model_v6.pdb
```

PyMOL ray trace modes:
- Mode 0: Default (clean, no shadows)
- Mode 1: Shadows (realistic, publication quality)
- Mode 2: Black outline only (line drawing)
- Mode 3: Quantized/posterized colors

## Common Data Variables

After sourcing `data_source_DE.R`, these key variables are available:
- `OGlcNAc_protein_DE_HEK293T/HepG2/Jurkat` - Differential expression results with logFC, adj.P.Val
- `OGlcNAc_protein_norm_HEK293T/HepG2/Jurkat` - Normalized TMT intensities
- `colors_cell`, `colors_glycan` - Consistent color palettes for plotting

## Analysis Logic Document

`OGlycoTM_Data_Analysis_Logic.RMD` contains the master analysis workflow organized by figure (Figure 1-6) with scheduled milestones.

## Key Data Insights

- **O-GlcNAc sites**: ~90% are in IDR (intrinsically disordered regions), ~10% in structured regions
- **IDR classification**: pLDDT < 50 = IDR, pLDDT ≥ 50 = Structured
- **Site features file**: `site_features/OGlcNAc_site_features.csv` contains per-site logFC, pLDDT, is_IDR, secondary structure
- **Proteins with sites in both regions**: Very rare (only 3-4 proteins across all cell types)
