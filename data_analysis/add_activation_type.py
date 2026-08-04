#!/usr/bin/env python3
"""
Add HCD/EThcD activation type to ranked O-GlcNAc PSM files.
This script reads calibrated mzML files to extract the activation method
(HCD or EThcD) for each spectrum and adds it to the ranked CSV files.

Uses indexed mzML access for fast direct scan lookup.
"""

import os
import pandas as pd
from pyteomics import mzml
from collections import defaultdict

# Configuration
SOURCE_PATH = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/"
OUTPUT_PATH = os.path.join(SOURCE_PATH, "point_to_point_response/")

# mzML directories for each cell type
MZML_DIRS = {
    "HEK293T": "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_HEK293T/",
    "HepG2": "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_HepG2/",
    "Jurkat": "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_Jurkat/"
}

# Ranked file paths
RANKED_FILES = {
    "HEK293T": os.path.join(OUTPUT_PATH, "OGlcNAc_Level1_HEK293T_ranked.csv"),
    "HepG2": os.path.join(OUTPUT_PATH, "OGlcNAc_Level1_HepG2_ranked.csv"),
    "Jurkat": os.path.join(OUTPUT_PATH, "OGlcNAc_Level1_Jurkat_ranked.csv")
}


def parse_spectrum_id(spectrum_str):
    """
    Parse spectrum identifier to extract file name and scan number.
    Format: {file_name}.{scan}.{scan}.{charge}
    Example: Eclipse_LF_OGlycoTM_HEK293T_OG_6_10062025_1uL.18782.18782.3

    Returns: (file_name, scan_number, charge)
    """
    parts = spectrum_str.rsplit('.', 3)
    if len(parts) == 4:
        file_name = parts[0]
        scan_number = int(parts[1])
        charge = int(parts[3])
        return file_name, scan_number, charge
    return None, None, None


def find_calibrated_mzml(file_name, mzml_dir):
    """
    Find the calibrated mzML file path for a given file name.
    Tries multiple naming conventions for calibrated files.
    """
    # Try different calibrated file naming patterns
    patterns = [
        f"{file_name}_calibrated.mzML",
        f"{file_name}_mz_calibrated.mzML",
        f"{file_name}_ppm_calibrated.mzML",
    ]

    for pattern in patterns:
        path = os.path.join(mzml_dir, pattern)
        if os.path.exists(path):
            return path

    # List all calibrated files and try to find a match
    try:
        all_files = os.listdir(mzml_dir)
        for f in all_files:
            if f.startswith(file_name) and 'calibrated' in f.lower() and f.endswith('.mzML'):
                return os.path.join(mzml_dir, f)
    except:
        pass

    return None


def get_activation_type_from_spectrum(spectrum):
    """
    Extract activation type from a spectrum dictionary.

    Returns: 'HCD', 'EThcD', 'ETD', or 'Unknown'
    """
    # Method 1: Check filter string (most reliable)
    if 'scanList' in spectrum and 'scan' in spectrum['scanList']:
        scan_info = spectrum['scanList']['scan'][0]
        filter_string = scan_info.get('filter string', '').lower()

        # EThcD has both @etd and @hcd in filter string
        # e.g., "873.1102@etd30.00 873.1102@hcd35.00"
        if '@etd' in filter_string and '@hcd' in filter_string:
            return 'EThcD'
        elif '@ethcd' in filter_string:
            return 'EThcD'
        elif '@etd' in filter_string:
            return 'ETD'
        elif '@hcd' in filter_string:
            return 'HCD'
        elif '@cid' in filter_string:
            return 'CID'

    # Method 2: Check activation dictionary in precursorList (fallback)
    if 'precursorList' in spectrum:
        precursor = spectrum['precursorList']['precursor'][0]
        if 'activation' in precursor:
            activation = precursor['activation']

            # Check for specific activation types
            if 'electron transfer dissociation' in activation:
                # Could be ETD or EThcD - need to check if HCD supplemental
                if 'beam-type collision-induced dissociation' in activation:
                    return 'EThcD'  # ETD with HCD supplemental
                return 'ETD'
            elif 'beam-type collision-induced dissociation' in activation:
                return 'HCD'
            elif 'collision-induced dissociation' in activation:
                return 'CID'

    return 'Unknown'


def get_activation_for_scan(mzml_reader, scan_number):
    """
    Get activation type for a specific scan using indexed access.

    Returns: activation type string or 'Unknown'
    """
    # Build the spectrum ID format used in mzML
    # Format: controllerType=0 controllerNumber=1 scan=XXXXX
    scan_id = f"controllerType=0 controllerNumber=1 scan={scan_number}"

    try:
        spectrum = mzml_reader.get_by_id(scan_id)
        return get_activation_type_from_spectrum(spectrum)
    except KeyError:
        # Try alternative ID format
        try:
            # Some mzML files use different ID formats
            for spec_id in mzml_reader.index.keys():
                if f"scan={scan_number}" in spec_id:
                    spectrum = mzml_reader.get_by_id(spec_id)
                    return get_activation_type_from_spectrum(spectrum)
        except:
            pass
    except Exception as e:
        pass

    return 'Unknown'


def add_activation_type_to_ranked_file(cell_type):
    """
    Add activation type column to a ranked CSV file for a given cell type.
    Uses indexed mzML access for fast lookup.
    """
    print(f"\n=== Processing {cell_type} ===")

    # Load ranked file
    ranked_file = RANKED_FILES[cell_type]
    if not os.path.exists(ranked_file):
        print(f"Ranked file not found: {ranked_file}")
        return None

    df = pd.read_csv(ranked_file)
    print(f"Loaded {len(df)} PSMs from ranked file")

    mzml_dir = MZML_DIRS[cell_type]

    # Group PSMs by mzML file to minimize file opens
    file_groups = defaultdict(list)
    for idx, row in df.iterrows():
        spectrum_str = row['Spectrum']
        file_name, scan_number, _ = parse_spectrum_id(spectrum_str)
        if file_name and scan_number:
            file_groups[file_name].append((idx, scan_number))

    print(f"Found PSMs from {len(file_groups)} unique mzML files")

    # Initialize activation type column
    df['Activation_Type'] = 'Unknown'

    # Process each mzML file with indexed access
    processed_files = 0
    total_scans_found = 0

    for file_name, scans in file_groups.items():
        mzml_path = find_calibrated_mzml(file_name, mzml_dir)

        if mzml_path is None:
            print(f"  Calibrated mzML not found for: {file_name}")
            continue

        print(f"  Opening: {os.path.basename(mzml_path)} ({len(scans)} scans)...", end=" ")

        try:
            # Use indexed mzML reader for fast random access
            with mzml.MzML(mzml_path, use_index=True) as reader:
                scans_found = 0
                for idx, scan_number in scans:
                    activation = get_activation_for_scan(reader, scan_number)
                    df.at[idx, 'Activation_Type'] = activation
                    if activation != 'Unknown':
                        scans_found += 1

                print(f"found {scans_found}/{len(scans)}")
                total_scans_found += scans_found
                processed_files += 1

        except Exception as e:
            print(f"Error: {e}")

    print(f"\nProcessed {processed_files}/{len(file_groups)} mzML files")
    print(f"Found activation type for {total_scans_found}/{len(df)} PSMs")

    # Summary of activation types
    print("\nActivation type summary:")
    activation_counts = df['Activation_Type'].value_counts()
    for act_type, count in activation_counts.items():
        print(f"  {act_type}: {count}")

    # Save updated file
    df.to_csv(ranked_file, index=False)
    print(f"\nSaved updated file: {ranked_file}")

    return df


if __name__ == "__main__":
    import sys

    # Process all cell types by default, or a specific one if provided
    if len(sys.argv) > 1:
        cell_types = [sys.argv[1]]
    else:
        cell_types = ["HEK293T", "HepG2", "Jurkat"]

    print("Adding activation type information to ranked O-GlcNAc PSM files")
    print("Using indexed mzML access for fast lookup")
    print("=" * 60)

    for cell_type in cell_types:
        if cell_type not in MZML_DIRS:
            print(f"Unknown cell type: {cell_type}")
            continue

        add_activation_type_to_ranked_file(cell_type)

    print("\n" + "=" * 60)
    print("Done!")
