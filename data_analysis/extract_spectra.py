#!/usr/bin/env python3
"""
Extract spectrum data from mzML files for top-ranked O-GlcNAc PSMs.
This script reads the ranked CSV files and extracts m/z and intensity data
from corresponding mzML files for spectral visualization.
"""

import os
import re
import pandas as pd
import numpy as np
from pyteomics import mzml

# Configuration
SOURCE_PATH = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/"
OUTPUT_PATH = os.path.join(SOURCE_PATH, "point_to_point_response/")

# mzML directories for each cell type
MZML_DIRS = {
    "HEK293T": "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_HEK293T/",
    "HepG2": "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_HepG2/",
    "Jurkat": "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_Jurkat/"
}


def parse_spectrum_id(spectrum_str):
    """
    Parse spectrum identifier to extract file name and scan number.
    Format: {file_name}.{scan}.{scan}.{charge}
    Example: Eclipse_LF_OGlycoTM_HEK293T_OG_6_10062025_1uL.18782.18782.3

    Returns: (file_name, scan_number, charge)
    """
    # Split by dots from the end (to handle file names with dots)
    parts = spectrum_str.rsplit('.', 3)
    if len(parts) == 4:
        file_name = parts[0]
        scan_number = int(parts[1])
        charge = int(parts[3])
        return file_name, scan_number, charge
    return None, None, None


def find_mzml_file(file_name, mzml_dir, use_calibrated=True):
    """
    Find the mzML file path for a given file name.
    Prefers calibrated version if use_calibrated=True.
    """
    if use_calibrated:
        calibrated_path = os.path.join(mzml_dir, f"{file_name}_calibrated.mzML")
        if os.path.exists(calibrated_path):
            return calibrated_path

    regular_path = os.path.join(mzml_dir, f"{file_name}.mzML")
    if os.path.exists(regular_path):
        return regular_path

    return None


def extract_spectrum(mzml_path, scan_number):
    """
    Extract spectrum data from mzML file for a specific scan number.

    Returns dict with:
        - mz_array: numpy array of m/z values
        - intensity_array: numpy array of intensity values
        - precursor_mz: precursor m/z
        - precursor_charge: precursor charge
        - ms_level: MS level (1 or 2)
        - retention_time: retention time in seconds
    """
    with mzml.read(mzml_path) as reader:
        for spectrum in reader:
            spec_id = spectrum.get('id', '')
            if f'scan={scan_number}' in spec_id:
                result = {
                    'mz_array': spectrum.get('m/z array', np.array([])),
                    'intensity_array': spectrum.get('intensity array', np.array([])),
                    'ms_level': spectrum.get('ms level'),
                    'spectrum_id': spec_id
                }

                # Get retention time from scanList
                if 'scanList' in spectrum and 'scan' in spectrum['scanList']:
                    scan_info = spectrum['scanList']['scan'][0]
                    result['retention_time'] = scan_info.get('scan start time', None)

                # Get precursor info for MS2 spectra
                if 'precursorList' in spectrum:
                    precursor = spectrum['precursorList']['precursor'][0]
                    if 'selectedIonList' in precursor:
                        sel_ion = precursor['selectedIonList']['selectedIon'][0]
                        result['precursor_mz'] = sel_ion.get('selected ion m/z')
                        result['precursor_charge'] = sel_ion.get('charge state')
                        result['precursor_intensity'] = sel_ion.get('peak intensity')

                return result

    return None


def extract_top_spectra(cell_type, n_top=50, use_calibrated=True):
    """
    Extract spectrum data for top N ranked PSMs from a cell type.

    Returns: DataFrame with spectrum info and file paths to extracted data
    """
    print(f"\n--- Extracting spectra for {cell_type} ---")

    # Load ranked file
    ranked_file = os.path.join(OUTPUT_PATH, f"OGlcNAc_Level1_{cell_type}_ranked.csv")
    df = pd.read_csv(ranked_file)
    print(f"Loaded {len(df)} PSMs")

    # Take top N
    top_df = df.head(n_top).copy()
    print(f"Processing top {len(top_df)} spectra")

    mzml_dir = MZML_DIRS[cell_type]

    # Create output directory for extracted spectra
    spectra_output_dir = os.path.join(OUTPUT_PATH, f"extracted_spectra_{cell_type}")
    os.makedirs(spectra_output_dir, exist_ok=True)

    results = []

    for idx, row in top_df.iterrows():
        spectrum_str = row['Spectrum']
        file_name, scan_number, charge = parse_spectrum_id(spectrum_str)

        if file_name is None:
            print(f"  Could not parse: {spectrum_str}")
            continue

        # Find mzML file
        mzml_path = find_mzml_file(file_name, mzml_dir, use_calibrated)

        if mzml_path is None:
            print(f"  mzML not found for: {file_name}")
            continue

        # Extract spectrum
        spec_data = extract_spectrum(mzml_path, scan_number)

        if spec_data is None:
            print(f"  Scan {scan_number} not found in {os.path.basename(mzml_path)}")
            continue

        # Save spectrum data to CSV
        spec_filename = f"{row['site_index']}_{scan_number}.csv"
        spec_filepath = os.path.join(spectra_output_dir, spec_filename)

        spec_df = pd.DataFrame({
            'mz': spec_data['mz_array'],
            'intensity': spec_data['intensity_array']
        })
        spec_df.to_csv(spec_filepath, index=False)

        # Collect result info
        results.append({
            'Quality_Rank': row['Quality_Rank'],
            'Gene': row['Gene'],
            'site_index': row['site_index'],
            'Peptide': row['Peptide'],
            'Composite_Score': row['Composite_Score'],
            'scan_number': scan_number,
            'precursor_mz': spec_data.get('precursor_mz'),
            'precursor_charge': spec_data.get('precursor_charge'),
            'ms_level': spec_data.get('ms_level'),
            'n_peaks': len(spec_data['mz_array']),
            'mzml_file': os.path.basename(mzml_path),
            'spectrum_file': spec_filename
        })

        if len(results) % 10 == 0:
            print(f"  Extracted {len(results)} spectra...")

    # Save results summary
    results_df = pd.DataFrame(results)
    summary_file = os.path.join(OUTPUT_PATH, f"extracted_spectra_{cell_type}_summary.csv")
    results_df.to_csv(summary_file, index=False)

    print(f"Extracted {len(results)} spectra")
    print(f"Summary saved to: {summary_file}")
    print(f"Spectra saved to: {spectra_output_dir}/")

    return results_df


if __name__ == "__main__":
    import sys

    # Default: extract top 50 from HEK293T
    cell_type = sys.argv[1] if len(sys.argv) > 1 else "HEK293T"
    n_top = int(sys.argv[2]) if len(sys.argv) > 2 else 50

    print(f"Extracting top {n_top} spectra from {cell_type}")
    results = extract_top_spectra(cell_type, n_top=n_top, use_calibrated=True)

    print("\n--- Top 10 extracted spectra ---")
    print(results[['Quality_Rank', 'Gene', 'site_index', 'n_peaks', 'precursor_mz']].head(10))
