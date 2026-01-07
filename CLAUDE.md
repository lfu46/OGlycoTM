# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

OGlycoTM is an R-based data analysis project for analyzing O-glycosylation (O-GlcNAc and O-GalNAc modifications) in mass spectrometry proteomics data across three human cell types: HEK293T, HepG2, and Jurkat. This is a research manuscript project generating publication-quality figures and statistical analysis.

## Data Pipeline Architecture

The analysis follows a sequential pipeline:

1. **data_import.R** - Imports raw TSV files from OGlycoTM mass spectrometry search results, filters for human proteins, removes decoys/contaminants
2. **data_filtering.R** - Filters for high-confidence localized sites (Level1, Level1b), removes cysteine artifacts, creates "bonafide" datasets
3. **data_source.R** - Central configuration with file paths, color palettes, and imports all filtered datasets
4. **Figure2.R** (and subsequent figure scripts) - Generates publication figures

## Running Scripts

```r
# Execute in order:
source('data_analysis/data_source.R')      # Load data and config (sources filtered data)
source('data_analysis/Figure2.R')          # Generate Figure 2 panels

# For full pipeline from raw data:
source('data_analysis/data_import.R')      # Import raw data
source('data_analysis/data_filtering.R')   # Filter & curate
```

## Key Configuration (data_source.R)

- **source_file_path**: `/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/`
- **figure_file_path**: `/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/`
- **Color palettes**:
  - Glycan types: O-GlcNAc (#F39B7F salmon), O-GalNAc (#4DBBD5 blue)
  - Cell types: HEK293T (#4DBBD5), HepG2 (#F39B7F), Jurkat (#00A087)

## Key R Packages

- `tidyverse` - Data manipulation and ggplot2 visualization
- `eulerr` - Proportional Venn/Euler diagrams
- `readr` - TSV/CSV I/O

## Data Files

All data files (CSV, TSV, XLSX) are in .gitignore. Data is stored on external network drive at `/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/`.

## Analysis Logic Document

`OGlycoTM_Data_Analysis_Logic.RMD` contains the master analysis workflow organized by figure (Figure 1-6) with scheduled milestones.
