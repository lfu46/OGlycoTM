#!/usr/bin/env python3
"""
Test script to regenerate example spectra with localization status.

Updates:
1. Shows "Localized" or "Not Localized" label
2. Correct precursor m/z from original data
3. Label positioned in upper right corner (not overlapping)
"""

import os
import sys
import numpy as np
import pandas as pd
from pyteomics import mzml
import matplotlib.pyplot as plt

# Add GlycoSpectrumAnnotator to path
sys.path.insert(0, '/Users/longpingfu/Downloads/GlycoSpectrumAnnotator')

from spectrum_annotator_ddzby import (
    SpectrumAnnotator,
    FragmentCalculator,
    parse_modifications_from_string,
)

# =============================================================================
# Configuration
# =============================================================================

SOURCE_PATH = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/"
OUTPUT_PATH = os.path.join(SOURCE_PATH, "low_prob_spectra/example_spectra/")

# Example spectra to regenerate
EXAMPLES = [
    {
        "name": "Example1_RBPJ_S38_unique",
        "cell_type": "HEK293T",
        "gene": "RBPJ",
        "peptide": "ANSSQVPPNESNTNSEGSYTNASTNSTSVTSSTATVVS",
        "is_localized": True,
        "description": "Localized"
    },
    {
        "name": "Example2_NUFIP2_S28_unique",
        "cell_type": "HEK293T",
        "gene": "NUFIP2",
        "peptide": "TIQNSSVSPTSSSSSSSSTGETQTQSSSR",
        "is_localized": True,
        "description": "Localized"
    },
    {
        "name": "Example3_ZNF362_S31_unique",
        "cell_type": "Jurkat",
        "gene": "ZNF362",
        "peptide": "TPSVSTSESSAGAGTGTGTSTPSTPTTTSQSR",
        "is_localized": True,
        "description": "Localized"
    },
    {
        "name": "Example4_UBN1_S10_weak",
        "cell_type": "Jurkat",
        "gene": "UBN1",
        "peptide": "LLTSPSLKPSAVSSVTSSTSLSK",
        "is_localized": False,
        "description": "Not Localized"
    },
]

# Bonafide files for getting precursor m/z
BONAFIDE_FILES = {
    "HEK293T": os.path.join(SOURCE_PATH, "filtered/OGlyco_HEK293T_bonafide.csv"),
    "HepG2": os.path.join(SOURCE_PATH, "filtered/OGlyco_HepG2_bonafide.csv"),
    "Jurkat": os.path.join(SOURCE_PATH, "filtered/OGlyco_Jurkat_bonafide.csv"),
}

# mzML directories
MZML_DIRS = {
    "HEK293T": "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_HEK293T/",
    "HepG2": "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_HepG2/",
    "Jurkat": "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_Jurkat/",
}


def find_calibrated_mzml(file_name, mzml_dir):
    """Find the calibrated mzML file path."""
    try:
        all_files = os.listdir(mzml_dir)
        for f in all_files:
            if f.startswith(file_name) and 'calibrated' in f.lower() and f.endswith('.mzML'):
                return os.path.join(mzml_dir, f)
    except:
        pass
    return None


def parse_spectrum_id(spectrum_str):
    """Parse spectrum identifier."""
    parts = spectrum_str.rsplit('.', 3)
    if len(parts) == 4:
        return parts[0], int(parts[1]), int(parts[3])
    return None, None, None


def extract_spectrum_data(mzml_reader, scan_number):
    """Extract spectrum data for a specific scan."""
    scan_id = f"controllerType=0 controllerNumber=1 scan={scan_number}"

    try:
        spectrum = mzml_reader.get_by_id(scan_id)
    except KeyError:
        for spec_id in mzml_reader.index.keys():
            if f"scan={scan_number}" in spec_id:
                spectrum = mzml_reader.get_by_id(spec_id)
                break
        else:
            return None

    return {
        'mz_array': spectrum.get('m/z array', np.array([])),
        'intensity_array': spectrum.get('intensity array', np.array([])),
    }


def add_localization_label(fig, is_localized):
    """Add localization status label to the figure."""
    if is_localized:
        label = "Localized"
        color = '#2E7D32'  # Green
        bg_color = '#E8F5E9'
    else:
        label = "Not Localized"
        color = '#C62828'  # Red
        bg_color = '#FFEBEE'

    # Create a fancy box
    bbox_props = dict(
        boxstyle="round,pad=0.3",
        facecolor=bg_color,
        edgecolor=color,
        linewidth=2
    )

    # Position at figure level - upper right corner (x=0.98, y=0.98 in figure coords)
    fig.text(0.98, 0.98, label,
             fontsize=10, fontweight='bold', fontfamily='Arial',
             color=color, verticalalignment='top', horizontalalignment='right',
             bbox=bbox_props, zorder=100)


def main():
    print("=" * 70)
    print("Regenerating Example Spectra with Localization Labels")
    print("=" * 70)

    # Load the annotation results
    results_file = os.path.join(SOURCE_PATH, "low_prob_spectra/low_prob_annotation_results.csv")
    if not os.path.exists(results_file):
        print(f"\nError: Results file not found: {results_file}")
        return

    df = pd.read_csv(results_file)
    print(f"\nLoaded {len(df)} annotation results")

    # Load bonafide files for precursor m/z
    bonafide_dfs = {}
    for cell_type, path in BONAFIDE_FILES.items():
        if os.path.exists(path):
            bonafide_dfs[cell_type] = pd.read_csv(path)
            print(f"  Loaded bonafide file for {cell_type}")

    os.makedirs(OUTPUT_PATH, exist_ok=True)

    # Process each example
    for example in EXAMPLES:
        print(f"\n{'='*70}")
        print(f"Processing: {example['name']}")
        print(f"  Gene: {example['gene']}")
        print(f"  Peptide: {example['peptide']}")
        print(f"  Status: {example['description']}")

        # Find matching PSM in results
        matches = df[
            (df['Gene'] == example['gene']) &
            (df['Peptide'] == example['peptide'])
        ]

        if len(matches) == 0:
            print(f"  WARNING: No matching PSM found!")
            continue

        # Use first match
        row = matches.iloc[0]
        print(f"  Found: {row['site_index']}, scan {row['scan_number']}")
        print(f"  Cell type: {row['Cell_Type']}")

        # Parse spectrum info
        spectrum_str = row['Spectrum']
        file_name, scan_number, charge = parse_spectrum_id(spectrum_str)

        # Get precursor m/z from bonafide file
        precursor_mz = 0
        if row['Cell_Type'] in bonafide_dfs:
            bf_match = bonafide_dfs[row['Cell_Type']][
                bonafide_dfs[row['Cell_Type']]['Spectrum'] == spectrum_str
            ]
            if len(bf_match) > 0:
                precursor_mz = bf_match.iloc[0].get('Observed.M.Z', 0)
                if pd.isna(precursor_mz):
                    precursor_mz = 0

        print(f"  Precursor m/z: {precursor_mz:.4f}")

        # Find mzML file
        mzml_dir = MZML_DIRS[row['Cell_Type']]
        mzml_path = find_calibrated_mzml(file_name, mzml_dir)

        if mzml_path is None:
            print(f"  WARNING: mzML file not found for {file_name}")
            continue

        print(f"  mzML: {os.path.basename(mzml_path)}")

        # Extract spectrum
        with mzml.MzML(mzml_path, use_index=True) as reader:
            spec_data = extract_spectrum_data(reader, scan_number)

        if spec_data is None:
            print(f"  WARNING: Could not extract spectrum")
            continue

        exp_mz = spec_data['mz_array']
        exp_intensity = spec_data['intensity_array']

        # Get modifications from bonafide file
        mod_string = ''
        if row['Cell_Type'] in bonafide_dfs:
            bf_match = bonafide_dfs[row['Cell_Type']][
                bonafide_dfs[row['Cell_Type']]['Spectrum'] == spectrum_str
            ]
            if len(bf_match) > 0:
                mod_string = bf_match.iloc[0].get('Assigned.Modifications', '')

        modifications = parse_modifications_from_string(mod_string if pd.notna(mod_string) else '')
        print(f"  Modifications: {mod_string}")

        # Find glycan position
        glycan_pos = None
        for mod in modifications:
            if abs(mod['mass'] - 528.2859) < 0.1:
                glycan_pos = mod['position']
                break

        print(f"  Glycan position: {glycan_pos}")

        # Create annotator
        try:
            annotator = SpectrumAnnotator(
                peptide=example['peptide'],
                modifications=modifications,
                precursor_charge=int(row['Charge']),
                precursor_mz=float(precursor_mz),
                exp_mz=exp_mz,
                exp_intensity=exp_intensity,
                tolerance_ppm=20.0,
                site_index=row['site_index'],
                gene=example['gene'],
                activation_type="HCD"
            )

            # Print matched ions with glycan markers
            print(f"\n  Matched b/y ions:")
            for ion in annotator.matched_ions:
                if ion.ion_type in ['b', 'y']:
                    glycan_marker = " (contains glycan)" if ion.has_modification else ""
                    print(f"    {ion.annotation}: m/z {ion.exp_mz:.4f}{glycan_marker}")

            # Generate figure without saving
            fig = annotator.plot(output_path=None)

            # Add localization label
            add_localization_label(fig, example['is_localized'])

            # Save figure
            output_file = os.path.join(OUTPUT_PATH, f"{example['name']}.pdf")
            fig.savefig(output_file, format='pdf', dpi=300, bbox_inches='tight')
            plt.close(fig)

            print(f"\n  Saved: {output_file}")

        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()

    print(f"\n{'='*70}")
    print("Done!")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
