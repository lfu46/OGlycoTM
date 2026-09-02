# MS/MS Spectrum Annotation (Python)

> Moved out of `CLAUDE.md` on 2026-09-02 so it is loaded on demand rather than into every session.

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

## False Match Rate Calculation

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

## Annotated spectrum folders (latest versions):
- **Tyr O-GlcNAc v5**: `Figures/Tyr_OGlcNAc_spectra_v5/{site}/` (4 selected sites: PRDX6_Y89, DDX17_Y580, HPRT1_Y105, SON_Y270) — MS1+HCD+EThcD+EMF
- **OGalNAc selected v2**: `Figures/OGalNAc_selected_spectra_v2/{site}/` (7 Jurkat sites) — MS1+HCD+EThcD
- **OGalNAc dropped protein check**: `Figures/OGalNAc_dropped_protein_check/` (11 proteins with zero/missing TMT)
- **OGalNAc nonsecretory check**: `Figures/OGalNAc_nonsecretory_check/` (101 HEK non-secretory EThcD PSMs) — MS1+HCD+EThcD
- **OGalNAc Level1/1b**: `Figures/OGalNAc_Level1_spectra/{cell}/` (1064 spectra), EThcD in `EThcD_only/`
- **ER/Golgi/PM v2**: `Figures/OGlcNAc_ER_Golgi_PM_spectra_v2/{cell}/{Level}/` (332 spectra), N-sequon in `N_sequon/`
- **Response figures**: `Figures/Response_Figures/` — Fig R1A (ER/Golgi/PM sequon analysis), Fig R1B (N-GlcNAc peptide overlap Venn: 6/289 = 2.1%)

## Spectrum annotator updates (Mar 29 2026):
- **c−NH3 removed**: c_n − NH3 = b_n identical mass, no longer generated as neutral loss
- **Clean labels**: single highest-priority label per peak, no "/" concatenation for overlapping ions
- **Y-ladder opt-in**: `extended_y_series=False` default globally; set `True` for natural complex glycans
- These changes are in the installed GlycoSpectrumAnnotator at `/Users/longpingfu/Downloads/GlycoSpectrumAnnotator/`

## O-GalNAc data (EThcD TMT fix applied at both site and protein level):
- Site PSMs: `data_source/site/OGalNAc_site_{cell}.csv` (originals backed up as `*_original.csv`)
- Filtered: `data_source/filtered/OGalNAc_{cell}.csv` (backups: `*_before_ethcd_fix.csv`)
- Site quant/norm: `data_source/quantification/OGalNAc_site_quant_{cell}.csv`, `normalization/OGalNAc_site_norm_{cell}.csv`
- Site DE (limma): `data_source/differential_analysis/OGalNAc_site_DE_{cell}.csv`
- Protein DE (limma): `data_source/differential_analysis/OGalNAc_protein_DE_{cell}.csv`
- Fix scripts: `fix_ogalnac_site_quant.py` (site), `regenerate_ogalnac_spectra.py` + inline (protein-level)

## O-GalNAc secretory classification:
- UniProt features: `data_source/reference/uniprot_features_OGalNAc.json` (485 proteins, full API JSON)
- Parsed features: `data_source/reference/uniprot_positional_features_OGalNAc.tsv`
- Site annotations: `data_source/annotation/OGalNAc_site_nielsen_regions.csv` (Nielsen et al. 2022 framework)
- **Secretory** = UniProt signal peptide OR transmembrane helix (following Steentoft EMBO J 2013)
- 31% secretory, 69% non-secretory across all O-GalNAc proteins — non-secretory enriched in O-GlcNAc functions
- Enrichment results: `data_source/enrichment/OGalNAc_{HepG2,Jurkat}_exclusive_GO.csv`, `OGalNAc_exclusive_GSEA_*_Jurkat.csv`
