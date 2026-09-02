# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

OGlycoTM is an R-based data analysis project for analyzing O-glycosylation (O-GlcNAc and O-GalNAc modifications) in mass spectrometry proteomics data across three human cell types: HEK293T, HepG2, and Jurkat. This is a research manuscript project generating publication-quality figures and statistical analysis.


## Reference material -- load on demand

Topic files in `docs/` (beside the static results site) are read when a task needs them rather than loaded into every session. Moved out of this file on 2026-09-02.

- Storage: `00_FILE_MAP.md` at the network project root is read before touching the data; artifact regeneration and the drive-drift checker -> `docs/storage_layout.md`
- PyMOL site panels (Figure 5E/5F), table-driven by `pymol_site_panels.py` + `.csv`, so a new panel is a CSV row and not a script -> `docs/pymol_structure_figures.md`
- MS/MS spectrum annotation with the installed `spectrum_annotator_ddzby` package: rules, FMR, spectrum folders, O-GalNAc data paths, secretory classification -> `docs/spectrum_annotation.md`
- Key data insights (IDR classification, site features file, OGlcNAc Atlas) -> `docs/key_data_insights.md`
- Document conversion recipes (Markdown/DOCX to PDF, PDF to TIFF, PDF to text) -> `docs/document_conversion.md`
- Revision and presentation notes, which point at the gitignored `REVISION_NOTES.local.md` to read before manuscript, figure-choice or reviewer questions -> `docs/revision_notes.md`

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
- **Env overrides**: set `OGLYCOTM_DATA` / `OGLYCOTM_FIGURES` to read from the faster Expansion
  working copy instead of the SMB share. Unset or unreadable falls back to the network path, so
  no script needs editing:
  ```r
  Sys.setenv(OGLYCOTM_DATA    = "/Volumes/Expansion/Longping/OGlycoTM/data_source")
  Sys.setenv(OGLYCOTM_FIGURES = "/Volumes/Expansion/Longping/OGlycoTM/Figures")
  ```
- **Color palettes**:
  - `colors_glycan`: O-GlcNAc (#F39B7F salmon), O-GalNAc (#4DBBD5 blue)
  - `colors_cell`: HEK293T (#4DBBD5), HepG2 (#F39B7F), Jurkat (#00A087)

## Storage layout (reorganized 2026-08-04)

The **network share is canonical** — 114 files hard-code `/Volumes/cos-lab-rwu60/...` with 232
references to the project root. Do not move or rename it.

| Where | Holds | Size |
|---|---|---|
| `/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/` | everything; `.raw` + calibrated `.mzML` live here only | ~161 GB |
| `/Volumes/Expansion/Longping/OGlycoTM/` | curated working copy: `data_source/`, `Figures/`, `Manuscript_Archive/`, repo snapshot | ~10 GB |
| this repo | code only | ~69 MB |


Internal revision notes are in `REVISION_NOTES.local.md`, gitignored because this repo is public.

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
- **Figure2.R** - O-GlcNAc identification overview (Euler diagrams, site distribution). **Output path: `Figure2_new/`** (not `Figure2/`). Has devEMF::emf() blocks for EMF export alongside PDF.
- **Figure3.R** - Differential expression analysis (volcano plots, heatmaps)
- **Figure4.R** - Cell-type specific analysis (circular heatmaps, GO enrichment). Has devEMF::emf() blocks for EMF export alongside PDF.
- **Figure5.R** - ~~Subcellular localization analysis~~ **REMOVED in Mar 2026 revision** (see Revision Notes below). Old Figure 6 becomes new Figure 5.
- **Figure6.R → Figure5.R (revised)** - Structural feature analysis (logFC distribution, secondary structure, IDR effects, site-specific examples)

**Figure delivery workflow (decided 2026-04-05):** Do NOT attempt EMF-for-editable-PPT. devEMF on macOS has multiple cross-platform issues (text as polygons, PLGBLT raster tiles for rotated elements, font metric mismatches). Deliver figures as TIFF/PDF. Change requests come back as text/markup; apply them by editing the R/Python source and re-exporting, never by editing the exported figure. Both Anal. Chem. and MCP accept TIFF (300 dpi RGB minimum), PDF, and EPS.

**CAUTION:** Do NOT `source("Figure2.R")` or `source("Figure4.R")` in full to regenerate a single panel — these scripts contain `enrichGO()` calls that take 5+ minutes. Instead, write a small driver script that loads the cached enrichment CSVs from `source_file_path/enrichment/` and rebuilds only the panel you need.

## Common Data Variables

After sourcing `data_source_DE.R`, these key variables are available:
- `OGlcNAc_protein_DE_HEK293T/HepG2/Jurkat` - Differential expression results with logFC, adj.P.Val
- `OGlcNAc_protein_norm_HEK293T/HepG2/Jurkat` - Normalized TMT intensities
- `colors_cell`, `colors_glycan` - Consistent color palettes for plotting

## Analysis Logic Document

`OGlycoTM_Data_Analysis_Logic.RMD` contains the master analysis workflow organized by figure (Figure 1-6) with scheduled milestones.
