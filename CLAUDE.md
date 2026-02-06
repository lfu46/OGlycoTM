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

## PyMOL 3D Structure Visualization (Figure 6E/6F)

PyMOL scripts for protein structure visualization with O-GlcNAc site highlighting:
- **Figure6E_*_pymol.py** - Candidate structured region sites (DDX50_Y492, PWP2_T23, PRDX6_T95, PRDX6_Y89)
- **Figure6F_*_pymol.py** - IDR region site examples (EWSR1_S274, etc.)

All PyMOL scripts share consistent formatting:
- pLDDT coloring (blue=high confidence, orange=low)
- Transparent surface (70% transparency)
- Site highlighted as colored spheres (orange=upregulated, cyan=stable)

```bash
# Run PyMOL script (requires PyMOL installed via Homebrew)
/opt/homebrew/bin/pymol -c -q Figure6E_PRDX6_T95_pymol.py

# AlphaFold structures are downloaded from:
# https://alphafold.ebi.ac.uk/files/AF-{UniProt_ID}-F1-model_v4.pdb (or v6)
```

PyMOL ray trace modes:
- Mode 0: Default (clean, no shadows)
- Mode 1: Shadows (realistic, publication quality)
- Mode 2: Black outline only (line drawing)
- Mode 3: Quantized/posterized colors

## MS/MS Spectrum Annotation (Python)

Python modules for EThcD spectrum extraction and annotation:
- **fragment_calculator.py** - Calculates theoretical m/z for b/y/c/z ions with modifications
- **spectrum_annotator.py** - Creates publication-quality annotated spectra (IPSA-style)
- **extract_ethcd_spectra.py** - Extracts EThcD spectra from calibrated mzML files

```bash
# Run spectrum annotation (requires pyteomics, numpy, pandas, matplotlib)
python spectrum_annotator.py

# Key parameters:
# - Mass tolerance: 20 ppm
# - Ion types: b, y (HCD), c, z (ETD), Y (intact glycopeptide), oxonium
# - Uses calibrated.mzML files for accurate mass matching
```

### False Match Rate Calculation

The spectrum annotator includes a false match rate (FMR) calculation based on the spectrum shifting method from Schulte et al. (Anal. Chem. 2025). This estimates the fraction of spurious matches:

```python
from fragment_calculator import calculate_false_match_rate

# Calculate FMR for a spectrum
fmr = calculate_false_match_rate(
    theoretical_ions,   # List of TheoreticalIon objects
    exp_mz,            # Experimental m/z array
    exp_intensity,     # Experimental intensity array
    tolerance_ppm=20.0,
    shift_range=25.0,  # Shift spectrum by π ± 25 Th
    shift_step=1.0     # 1 Th increments
)

print(f"FMR (peaks): {fmr.fmr_peaks*100:.1f}%")
print(f"FMR (intensity): {fmr.fmr_intensity*100:.1f}%")
```

The method shifts the spectrum by π ± 25 Th (π offset prevents isotope pattern matches) and calculates what fraction of matches would occur by chance. Low FMR (<10%) indicates high-quality annotations.

Spectrum data locations:
- mzML files: `/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_{cell_type}/`
- EThcD ranked files: `data_source/point_to_point_response/OGlcNAc_Level1_{cell_type}_EThcD_ranked.csv`
- Extracted spectra: `data_source/point_to_point_response/{cell_type}_ethcd_spectra/`

## Common Data Variables

After sourcing `data_source_DE.R`, these key variables are available:
- `OGlcNAc_protein_DE_HEK293T/HepG2/Jurkat` - Differential expression results with logFC, adj.P.Val
- `OGlcNAc_protein_norm_HEK293T/HepG2/Jurkat` - Normalized TMT intensities
- `colors_cell`, `colors_glycan` - Consistent color palettes for plotting

## Analysis Logic Document

`OGlycoTM_Data_Analysis_Logic.RMD` contains the master analysis workflow organized by figure (Figure 1-6) with scheduled milestones.

## Key Data Insights

- **O-GlcNAc sites**: ~90% are in IDR (intrinsically disordered regions), ~10% in structured regions
- **IDR classification**: Uses StructureMap's pPSE method - smoothed `nAA_24_180_pae` (pPSE_24_smooth10) ≤ 34.27 = IDR. This is calculated in `structuremap_analysis.py` following the StructureMap tutorial. Note: pLDDT is included as a regression predictor but is NOT used for IDR classification.
- **Site features file**: `site_features/OGlcNAc_site_features.csv` contains per-site logFC, pLDDT, is_IDR, secondary structure, pPSE_24, pPSE_12, pPSE_24_smooth10
- **Proteins with sites in both regions**: Very rare (only 3-4 proteins across all cell types)
- **OGlcNAc Atlas reference**: `reference/OGlcNAcAtlas_unambiguous_sites_20251208.csv` for checking reported vs novel sites

## Markdown to PDF Conversion

Convert markdown files to text-selectable PDF using pandoc:

```bash
pandoc "input.md" -o "output.pdf" --pdf-engine=xelatex \
  -V geometry:margin=1in -V fontsize=11pt -V mainfont="Times New Roman" \
  --toc -V colorlinks=true -V linkcolor=blue -V urlcolor=blue
```

If Unicode subscripts/superscripts cause issues (Times New Roman doesn't support them), pipe through sed first:
```bash
sed -e 's/⁻/-/g' -e 's/₀/0/g' -e 's/₁/1/g' -e 's/₂/2/g' -e 's/₃/3/g' -e 's/₄/4/g' -e 's/₅/5/g' -e 's/₆/6/g' -e 's/₇/7/g' -e 's/₈/8/g' -e 's/₉/9/g' "input.md" | \
  pandoc -o "output.pdf" --pdf-engine=xelatex \
  -V geometry:margin=1in -V fontsize=11pt -V mainfont="Times New Roman" \
  --toc -V colorlinks=true -V linkcolor=blue -V urlcolor=blue
```

## DOCX to PDF Conversion

Convert Word documents to text-selectable PDF using pandoc:

```bash
pandoc "input.docx" -o "output.pdf" --pdf-engine=xelatex
```

Note: This produces text-selectable PDFs but may not preserve complex formatting (tables, images, custom styles) perfectly. For exact formatting preservation, use `docx2pdf` (`pip install docx2pdf`) which requires Microsoft Word or LibreOffice.

## PDF to High-Resolution TIFF Conversion

Convert PDF files to high-resolution TIFF (600 DPI) using ImageMagick:

```bash
# Single file
magick -density 600 "input.pdf" -quality 100 "output.tiff"

# Batch convert all PDFs in a folder
for pdf in /path/to/folder/*.pdf; do
  filename=$(basename "$pdf" .pdf)
  magick -density 600 "$pdf" -quality 100 "/path/to/output/${filename}.tiff"
done
```

Note: Requires ImageMagick (`brew install imagemagick`). 600 DPI TIFFs are large (~80-90 MB each) but publication-quality.

## PDF to Text Extraction (Token-Efficient)

Extract text from PDF files for uploading to Claude Desktop (avoids image tokens):

```bash
# Requires poppler: brew install poppler

# Single file (compact, token-efficient)
pdftotext "input.pdf" "output.txt"

# Batch convert all PDFs in a folder
for pdf in /path/to/folder/*.pdf; do
  pdftotext "$pdf" "${pdf%.pdf}.txt"
done
```

Note: The default mode (no flags) is most token-efficient. Use `-layout` only if you need to preserve table/column structure (costs more tokens due to extra whitespace).
