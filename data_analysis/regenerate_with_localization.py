#!/usr/bin/env python3
"""
Regenerate all Level1b spectra with localization status annotations.

This script:
1. Loads the site localization analysis results
2. Regenerates each spectrum with "Localized" or "Not Localized" label
3. Combines all spectra into a single PDF
"""

import os
import sys
import json
import numpy as np
import pandas as pd
from pyteomics import mzml
from collections import defaultdict
import matplotlib.pyplot as plt
from PyPDF2 import PdfMerger

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
OUTPUT_PATH = os.path.join(SOURCE_PATH, "low_prob_spectra/")

# Input files
ANNOTATION_RESULTS = os.path.join(OUTPUT_PATH, "low_prob_annotation_results.csv")
LOCALIZATION_ANALYSIS = os.path.join(OUTPUT_PATH, "site_localization_analysis.csv")

# Bonafide files for modifications and precursor m/z
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


def parse_spectrum_id(spectrum_str):
    """Parse spectrum identifier."""
    parts = spectrum_str.rsplit('.', 3)
    if len(parts) == 4:
        return parts[0], int(parts[1]), int(parts[3])
    return None, None, None


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

    # Position at figure level - upper right corner
    fig.text(0.98, 0.98, label,
             fontsize=10, fontweight='bold', fontfamily='Arial',
             color=color, verticalalignment='top', horizontalalignment='right',
             bbox=bbox_props, zorder=100)


def main():
    print("=" * 70)
    print("Regenerating Spectra with Localization Annotations")
    print("=" * 70)

    # Load data
    print("\nLoading data...")
    annotation_df = pd.read_csv(ANNOTATION_RESULTS)
    localization_df = pd.read_csv(LOCALIZATION_ANALYSIS)

    print(f"  Annotation results: {len(annotation_df)} PSMs")
    print(f"  Localization analysis: {len(localization_df)} PSMs")

    # Create localization lookup by site_index + peptide
    loc_info = {}
    for _, row in localization_df.iterrows():
        key = row['site_index'] + '_' + row['Peptide']
        # Localized = unique or single_site
        is_localized = row['status'] in ['unique', 'single_site']
        loc_info[key] = is_localized

    # Load bonafide files for modifications and precursor m/z
    bonafide_dfs = {}
    for cell_type, path in BONAFIDE_FILES.items():
        if os.path.exists(path):
            bonafide_dfs[cell_type] = pd.read_csv(path)
            print(f"  Loaded bonafide file for {cell_type}")

    # Create output directory
    output_dir = os.path.join(OUTPUT_PATH, "annotated_spectra_with_status")
    os.makedirs(output_dir, exist_ok=True)

    # Group by cell type and file
    file_groups = defaultdict(lambda: defaultdict(list))
    for idx, row in annotation_df.iterrows():
        cell_type = row['Cell_Type']
        spectrum_str = row['Spectrum']
        file_name, scan_number, _ = parse_spectrum_id(spectrum_str)
        if file_name and scan_number:
            file_groups[cell_type][file_name].append((idx, scan_number, row))

    # Process each spectrum
    n_processed = 0
    n_localized = 0
    n_not_localized = 0
    n_errors = 0
    output_files = []

    for cell_type in ['HEK293T', 'HepG2', 'Jurkat']:
        if cell_type not in file_groups:
            continue

        mzml_dir = MZML_DIRS[cell_type]
        files = file_groups[cell_type]

        print(f"\n{cell_type}: {sum(len(scans) for scans in files.values())} PSMs")

        for file_name, scans in files.items():
            mzml_path = find_calibrated_mzml(file_name, mzml_dir)

            if mzml_path is None:
                print(f"  WARNING: mzML not found: {file_name}")
                n_errors += len(scans)
                continue

            print(f"  Processing: {os.path.basename(mzml_path)} ({len(scans)} scans)...", end=" ", flush=True)

            try:
                with mzml.MzML(mzml_path, use_index=True) as reader:
                    scans_processed = 0

                    for idx, scan_number, row in scans:
                        spec_data = extract_spectrum_data(reader, scan_number)

                        if spec_data is None:
                            n_errors += 1
                            continue

                        exp_mz = spec_data['mz_array']
                        exp_intensity = spec_data['intensity_array']

                        if len(exp_mz) == 0:
                            n_errors += 1
                            continue

                        # Get modifications and precursor m/z from bonafide file
                        mod_string = ''
                        precursor_mz = 0
                        spectrum_str = row['Spectrum']

                        if cell_type in bonafide_dfs:
                            bf_match = bonafide_dfs[cell_type][
                                bonafide_dfs[cell_type]['Spectrum'] == spectrum_str
                            ]
                            if len(bf_match) > 0:
                                mod_string = bf_match.iloc[0].get('Assigned.Modifications', '')
                                precursor_mz = bf_match.iloc[0].get('Observed.M.Z', 0)
                                if pd.isna(precursor_mz):
                                    precursor_mz = 0
                                if pd.isna(mod_string):
                                    mod_string = ''

                        modifications = parse_modifications_from_string(mod_string)

                        peptide = row['Peptide']
                        precursor_charge = int(row['Charge'])

                        # Get localization status
                        pep_key = row['site_index'] + '_' + peptide
                        is_localized = loc_info.get(pep_key, False)

                        if is_localized:
                            n_localized += 1
                        else:
                            n_not_localized += 1

                        try:
                            # Create annotator
                            annotator = SpectrumAnnotator(
                                peptide=peptide,
                                modifications=modifications,
                                precursor_charge=precursor_charge,
                                precursor_mz=float(precursor_mz),
                                exp_mz=exp_mz,
                                exp_intensity=exp_intensity,
                                tolerance_ppm=20.0,
                                site_index=row['site_index'],
                                gene=row['Gene'],
                                activation_type="HCD"
                            )

                            # Create figure without saving
                            fig = annotator.plot(output_path=None)

                            # Add localization label
                            add_localization_label(fig, is_localized)

                            # Save figure
                            output_file = os.path.join(output_dir, f"{row['site_index']}_{scan_number}.pdf")
                            fig.savefig(output_file, format='pdf', dpi=300, bbox_inches='tight')
                            plt.close(fig)

                            output_files.append(output_file)
                            scans_processed += 1
                            n_processed += 1

                        except Exception as e:
                            n_errors += 1
                            continue

                    print(f"done ({scans_processed})")

            except Exception as e:
                print(f"Error: {e}")
                n_errors += len(scans)

    print(f"\n{'='*70}")
    print(f"Processing Complete")
    print(f"{'='*70}")
    print(f"  Processed: {n_processed}")
    print(f"  Localized: {n_localized}")
    print(f"  Not Localized: {n_not_localized}")
    print(f"  Errors: {n_errors}")

    # Combine all PDFs
    print(f"\nCombining {len(output_files)} PDFs...")

    merger = PdfMerger()

    # Sort files for consistent ordering
    output_files_sorted = sorted(output_files)

    for i, pdf_file in enumerate(output_files_sorted):
        if (i + 1) % 100 == 0:
            print(f"  Adding file {i + 1}/{len(output_files_sorted)}...")
        try:
            merger.append(pdf_file)
        except Exception as e:
            print(f"  Error adding {os.path.basename(pdf_file)}: {e}")

    combined_output = os.path.join(OUTPUT_PATH, "all_annotated_spectra_with_status.pdf")
    merger.write(combined_output)
    merger.close()

    file_size = os.path.getsize(combined_output) / (1024 * 1024)
    print(f"\nCombined PDF saved to:")
    print(f"  {combined_output}")
    print(f"  Size: {file_size:.1f} MB")
    print(f"  Pages: {len(output_files_sorted)}")


if __name__ == "__main__":
    main()
