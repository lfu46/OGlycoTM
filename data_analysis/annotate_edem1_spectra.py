#!/usr/bin/env python3
"""Annotate all EDEM1 O-GalNAc EThcD spectra using GlycoSpectrumAnnotator."""

import numpy as np
import os
import sys

import mzml_utils
from spectrum_annotator_ddzby import SpectrumAnnotator
from spectrum_annotator_ddzby.fragment_calculator import parse_modifications_from_string

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Arial'

# Output directory
OUT_DIR = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/EDEM1_OGalNAc_spectra"
os.makedirs(OUT_DIR, exist_ok=True)

# Base path for mzML files
BASE = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version"

# All EDEM1 O-GalNAc PSMs
PSMS = [
    # (cell_type, mzml_basename, scan, charge, confidence, prob, mod_string)
    ("HEK293T", "Eclipse_LF_OGlycoTM_HEK293T_OG_1_10062025", 16906, 3, "Level1b", 0.328,
     "N-term(229.1629),1S(555.2968)"),
    ("HepG2", "Eclipse_LF_OGlycoTM_HepG2_OG_1_11212025", 24425, 3, "Level1", 0.943,
     "N-term(229.1629),1S(555.2968)"),
    ("HepG2", "Eclipse_LF_OGlycoTM_HepG2_OG_1_11212025", 24433, 3, "Level1", 0.916,
     "N-term(229.1629),1S(555.2968)"),
    ("Jurkat", "Eclipse_LF_OGlycoTM_Jurkat_OG_1_10222025", 19300, 3, "Level1", 0.884,
     "N-term(229.1629),1S(555.2968)"),
    ("Jurkat", "Eclipse_LF_OGlycoTM_Jurkat_OG_1_10222025", 19749, 3, "Level1", 0.944,
     "N-term(229.1629),1S(555.2968)"),
    ("Jurkat", "Eclipse_LF_OGlycoTM_Jurkat_OG_24_10222025", 17548, 3, "Level1", 0.875,
     "N-term(229.1629),1S(555.2968)"),
    ("Jurkat", "Eclipse_LF_OGlycoTM_Jurkat_OG_2_10222025", 18679, 3, "Level1", 0.944,
     "N-term(229.1629),1S(555.2968)"),
    ("Jurkat", "Eclipse_LF_OGlycoTM_Jurkat_OG_2_10222025", 18690, 3, "Level1", 0.940,
     "N-term(229.1629),1S(555.2968)"),
]

PEPTIDE = "SPDGPASPTSGPVGR"

# Cache opened readers
readers = {}

for cell, basename, scan, charge, conf, prob, mod_string in PSMS:
    print(f"\n=== {cell} scan {scan} ({conf}, prob={prob:.3f}) ===")

    # Find and open mzML
    mzml_dir = os.path.join(BASE, f"OGlycoTM_{cell}")
    mzml_path = os.path.join(mzml_dir, f"{basename}_calibrated.mzML")

    if mzml_path not in readers:
        print(f"  Opening {os.path.basename(mzml_path)}...")
        readers[mzml_path] = mzml_utils.MzMLReader(mzml_path)
    reader = readers[mzml_path]

    # Get spectrum
    spec = reader.get_spectrum(scan)
    if spec is None:
        print(f"  ERROR: scan {scan} not found")
        continue

    # Check activation type from filter string
    fs = spec.filter_string if hasattr(spec, 'filter_string') else ""
    if "etd" in fs.lower() or "ethcd" in fs.lower():
        activation = "EThcD"
    elif "hcd" in fs.lower():
        activation = "HCD"
    else:
        activation = "EThcD"  # default for OPair

    print(f"  Filter: {fs}")
    print(f"  Activation: {activation}")
    print(f"  Precursor m/z: {spec.precursor_mz:.4f}, peaks: {spec.n_peaks}")

    # Parse modifications
    mods = parse_modifications_from_string(mod_string)

    # Annotate
    annotator = SpectrumAnnotator(
        peptide=PEPTIDE,
        modifications=mods,
        precursor_charge=charge,
        precursor_mz=spec.precursor_mz,
        exp_mz=spec.mz,
        exp_intensity=spec.intensity,
        tolerance_ppm=20.0,
        activation_type=activation,
        do_deisotope=False,
        scan_num=scan,
        sn_threshold=0.0,
        confidence_level=conf,
        site_index="Q92611_S48",
        gene="EDEM1",
        source_file=basename,
    )

    fig = annotator.plot()

    # Save
    tag = f"EDEM1_S48_{cell}_scan{scan}_{activation}"
    pdf_path = os.path.join(OUT_DIR, f"{tag}.pdf")
    png_path = os.path.join(OUT_DIR, f"{tag}.png")
    fig.savefig(pdf_path, dpi=600, bbox_inches='tight')
    fig.savefig(png_path, dpi=600, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {tag}.pdf")

print(f"\nAll spectra saved to: {OUT_DIR}")
print(f"Files: {sorted(os.listdir(OUT_DIR))}")
