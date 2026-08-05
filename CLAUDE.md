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

Read **`00_FILE_MAP.md` at the network project root** before touching the data. It records the
layout, the search provenance (FragPipe 24.0 / MSFragger 4.4 / pwiz 3.0.25139), and the exact
commands that regenerate the two deleted classes of artifact:

- 43 GB of 600-DPI TIFF spectra — rebuilt from the kept PDFs with `magick -density 600 … -quality 100`
  (verified pixel-identical).
- 44.7 GB of **uncalibrated** `.mzML` — rebuilt from the kept `.raw` with the recorded `msconvert`
  filter chain. Analysis uses the *calibrated* mzML, which was kept.

Drift between the two drives is declared in `oglycotm_sync_policy.tsv` on the Expansion drive and
checked with the shared checker:

```bash
python3 /Users/longpingfu/Downloads/OGlyco_DBA/data_analysis/tools/dba_sync_check.py \
        --policy /Volumes/Expansion/Longping/OGlycoTM/oglycotm_sync_policy.tsv
#   ... --checksum   md5-verify same-size files (catches SMB truncation)
#   ... --apply      reconcile (never deletes; a size mismatch needs a human)
```

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

## PyMOL 3D Structure Visualization (Figure 5E/5F, formerly 6E/6F)

Site panels are **table-driven**. The 16 `OGalNAc_*_pymol.py` and 4 `Figure6E_*_pymol.py` scripts
were two templates repeated with three values changed; they were replaced 2026-08-04 by
`pymol_site_panels.py` + `pymol_site_panels.csv` (2,143 lines → 309), verified to render
pixel-identically.

```bash
python3 pymol_site_panels.py --list          # show the panel table
python3 pymol_site_panels.py --only PTPRC    # render one (matches gene/accession/site)
python3 pymol_site_panels.py                 # render all 20
```

To add a panel, add a CSV row — do not write a new script. Columns:
`style,gene,accession,sites,color,out_subdir,prefix,notes`, where `style` is `surface`
(whole residue as spheres, the O-GalNAc panels) or `sidechain` (only the modified side-chain
atoms, the Figure 6E panels), and `sites` is a semicolon list like `T139;S146`.

Shared formatting (unchanged from the originals):
- pLDDT coloring (blue=high confidence, orange=low)
- Transparent surface (70% transparency)
- Site highlighted as colored spheres (orange=upregulated, cyan=stable)

**NEVER hand-write an AlphaFold filename or version.** The driver resolves structures via
`mzml_utils.structure.fetch_structure()` with a cached-file fast path. Interpolating
`AF-{acc}-F1-model_v6.pdb` is how `Figure6F_EWSR1_S274_pymol.py` ended up pinned to `model_v4`
while every other panel used v6 — that script now resolves the newest cached model instead.

Still hand-written, kept for genuine per-protein domain/IDR colouring that does not parameterize:
`Figure6F_EWSR1_S274_pymol.py`, `Figure6F_HOXA13_pymol.py`, `Figure6F_HYOU1_pymol.py`,
`Figure6F_pymol.py`, and `Figure6F_ray_modes.py` (render-settings exploration).

```bash
# Run a hand-written PyMOL script (requires PyMOL installed via Homebrew)
/opt/homebrew/bin/pymol -c -q Figure6F_HYOU1_pymol.py
```

PyMOL ray trace modes:
- Mode 0: Default (clean, no shadows)
- Mode 1: Shadows (realistic, publication quality)
- Mode 2: Black outline only (line drawing)
- Mode 3: Quantized/posterized colors

## MS/MS Spectrum Annotation (Python)

**CRITICAL: ALWAYS read spectra via `mzml_utils.open_spectra(path)` — NEVER `pyteomics.mzml` (streaming, extremely slow), and do not construct `MzMLReader` directly. `open_spectra` is a drop-in that returns a fast `SpectrumCache` when `<mzml_dir>/spectra_cache/<stem>.spectra.db` exists and an `MzMLReader` otherwise, so it is never slower and gets faster the moment a cache is built. Every active script here was converted 2026-08-04. Always check for cached data (pickles, filter_string_cache.csv) before opening mzML files.**

**Use GlycoSpectrumAnnotator** (`spectrum_annotator_ddzby`, installed as editable package at `/Users/longpingfu/Downloads/GlycoSpectrumAnnotator/`) for all annotation. The local `spectrum_annotator.py` / `fragment_calculator.py` forks were **deleted 2026-08-04** — they were strict older subsets (646 vs 1490 and 1275 vs 1855 lines), missing the glycan library, N-glycan support, precursor-envelope filtering and the isotope-consistency flags. Every script behind a published spectrum already used the installed package. Import from `spectrum_annotator_ddzby` / `spectrum_annotator_ddzby.fragment_calculator`, never from a local module.

Python modules for spectrum annotation:
- **mzml_utils** (`import mzml_utils`) - Cache-aware spectrum reader (`open_spectra`), ion search, fragment calculator, deisotoping, spectral similarity, protease digestion
- **GlycoSpectrumAnnotator** (`spectrum_annotator_ddzby`) - Publication-quality annotated spectra with correct butterfly diagram, glycan labels, deisotoping, S/N filtering, charge-reduced exclusion
- **OGlyco_DBA tools** (`/Users/longpingfu/Downloads/OGlyco_DBA/data_analysis/`) - `opair_utils.py`, `oglyco_validation.py`, `mass_degeneracy.py` for validation workflows
- **extract_ethcd_spectra.py** - Extracts EThcD spectra from calibrated mzML files

```python
from spectrum_annotator_ddzby import SpectrumAnnotator
from spectrum_annotator_ddzby.fragment_calculator import load_noise_cache

# Key parameters for annotation:
annotator = SpectrumAnnotator(
    peptide=peptide, modifications=mods,
    precursor_charge=charge, precursor_mz=obs_mz,
    exp_mz=mz, exp_intensity=intensity,
    tolerance_ppm=20.0,
    activation_type='EThcD',  # or 'HCD'
    do_deisotope=False,       # NO MS2 deisotoping (matches MSFragger behavior)
    scan_num=scan,
    sn_threshold=0.0,         # 0.0 for manual validation, 5.0 for automated
    confidence_level='Level1', # shown in title
)
```

**Critical annotation rules (updated 2026-03-26):**
- Use **calibrated mzML** files for annotation (MSFragger output), NOT raw-converted mzML
- OPair modification masses: check `Assigned.Modifications` per PSM (299.123 vs 528.286 for HexNAc)
- CAM: only add if explicitly listed in `Assigned.Modifications`
- **Bare b/y ions cannot localize glycosites in HCD** — need glycan-retaining ions (EThcD c/z)
- For manual validation: MS1 (isolation + mass accuracy) → HCD (oxonium + sequence) → EThcD (localization)
- **z-ions use z+H (even-electron) formula**, NOT z• (radical). Matches MSFragger. Δ = 1.007 Da.
- **No MS2 deisotoping** (`do_deisotope=False`) — matches MSFragger behavior
- **No isotope matching** — prevents c/z complementary ion overlap artifacts
- **Residue-specific neutral losses**: H2O only for S/T/D/E; NH3 only for R/K/N/Q; CO2 removed
- **EThcD TMT quantification**: use paired HCD scan reporters (same precursor, higher NCE yield). Script: `fix_ogalnac_site_quant.py`

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

### Annotated spectrum folders (latest versions):
- **Tyr O-GlcNAc v5**: `Figures/Tyr_OGlcNAc_spectra_v5/{site}/` (4 selected sites: PRDX6_Y89, DDX17_Y580, HPRT1_Y105, SON_Y270) — MS1+HCD+EThcD+EMF
- **OGalNAc selected v2**: `Figures/OGalNAc_selected_spectra_v2/{site}/` (7 Jurkat sites) — MS1+HCD+EThcD
- **OGalNAc dropped protein check**: `Figures/OGalNAc_dropped_protein_check/` (11 proteins with zero/missing TMT)
- **OGalNAc nonsecretory check**: `Figures/OGalNAc_nonsecretory_check/` (101 HEK non-secretory EThcD PSMs) — MS1+HCD+EThcD
- **OGalNAc Level1/1b**: `Figures/OGalNAc_Level1_spectra/{cell}/` (1064 spectra), EThcD in `EThcD_only/`
- **ER/Golgi/PM v2**: `Figures/OGlcNAc_ER_Golgi_PM_spectra_v2/{cell}/{Level}/` (332 spectra), N-sequon in `N_sequon/`
- **Response figures**: `Figures/Response_Figures/` — Fig R1A (ER/Golgi/PM sequon analysis), Fig R1B (N-GlcNAc peptide overlap Venn: 6/289 = 2.1%)

### Spectrum annotator updates (Mar 29 2026):
- **c−NH3 removed**: c_n − NH3 = b_n identical mass, no longer generated as neutral loss
- **Clean labels**: single highest-priority label per peak, no "/" concatenation for overlapping ions
- **Y-ladder opt-in**: `extended_y_series=False` default globally; set `True` for natural complex glycans
- These changes are in the installed GlycoSpectrumAnnotator at `/Users/longpingfu/Downloads/GlycoSpectrumAnnotator/`

### O-GalNAc data (EThcD TMT fix applied at both site and protein level):
- Site PSMs: `data_source/site/OGalNAc_site_{cell}.csv` (originals backed up as `*_original.csv`)
- Filtered: `data_source/filtered/OGalNAc_{cell}.csv` (backups: `*_before_ethcd_fix.csv`)
- Site quant/norm: `data_source/quantification/OGalNAc_site_quant_{cell}.csv`, `normalization/OGalNAc_site_norm_{cell}.csv`
- Site DE (limma): `data_source/differential_analysis/OGalNAc_site_DE_{cell}.csv`
- Protein DE (limma): `data_source/differential_analysis/OGalNAc_protein_DE_{cell}.csv`
- Fix scripts: `fix_ogalnac_site_quant.py` (site), `regenerate_ogalnac_spectra.py` + inline (protein-level)

### O-GalNAc secretory classification:
- UniProt features: `data_source/reference/uniprot_features_OGalNAc.json` (485 proteins, full API JSON)
- Parsed features: `data_source/reference/uniprot_positional_features_OGalNAc.tsv`
- Site annotations: `data_source/annotation/OGalNAc_site_nielsen_regions.csv` (Nielsen et al. 2022 framework)
- **Secretory** = UniProt signal peptide OR transmembrane helix (following Steentoft EMBO J 2013)
- 31% secretory, 69% non-secretory across all O-GalNAc proteins — non-secretory enriched in O-GlcNAc functions
- Enrichment results: `data_source/enrichment/OGalNAc_{HepG2,Jurkat}_exclusive_GO.csv`, `OGalNAc_exclusive_GSEA_*_Jurkat.csv`

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

## Revision / presentation notes

The Anal. Chem. revision record and cohort-talk notes live in `REVISION_NOTES.local.md`, which is
gitignored because this repository is public. Read that file before answering manuscript,
figure-choice, or reviewer-response questions.

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
