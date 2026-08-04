#!/usr/bin/env python3
"""
Extract activation type (HCD vs ETD/EThcD) from mzML files for low probability PSMs.
Uses indexed mzML access for fast extraction.
"""

import os
import pandas as pd
from pyteomics import mzml
from collections import defaultdict, Counter

# =============================================================================
# Configuration
# =============================================================================

# Input: Low probability PSM file
LOW_PROB_PSM_FILE = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/site/OGlcNAc_low_prob_psm.csv"

# mzML directories (using calibrated files)
MZML_DIRS = {
    "HEK293T": "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_HEK293T/",
    "HepG2": "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_HepG2/",
    "Jurkat": "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_Jurkat/"
}

# Output file
OUTPUT_FILE = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/site/OGlcNAc_low_prob_psm_activation.csv"


# =============================================================================
# Helper functions
# =============================================================================

def parse_spectrum_id(spectrum_str):
    """
    Parse spectrum identifier to extract file name and scan number.
    Format: {file_name}.{scan}.{scan}.{charge}
    """
    parts = spectrum_str.rsplit('.', 3)
    if len(parts) == 4:
        file_name = parts[0]
        scan_number = int(parts[1])
        charge = int(parts[3])
        return file_name, scan_number, charge
    return None, None, None


def find_calibrated_mzml(file_name, mzml_dir):
    """Find the calibrated mzML file path for a given file name."""
    patterns = [
        f"{file_name}_calibrated.mzML",
        f"{file_name}_mz_calibrated.mzML",
        f"{file_name}_ppm_calibrated.mzML",
    ]

    for pattern in patterns:
        path = os.path.join(mzml_dir, pattern)
        if os.path.exists(path):
            return path

    # Search for any matching calibrated file
    try:
        all_files = os.listdir(mzml_dir)
        for f in all_files:
            if f.startswith(file_name) and 'calibrated' in f.lower() and f.endswith('.mzML'):
                return os.path.join(mzml_dir, f)
    except:
        pass

    return None


def get_activation_type(mzml_reader, scan_number):
    """
    Extract activation type for a specific scan using indexed access.
    Returns: 'HCD', 'ETD', 'EThcD', or 'Unknown'
    """
    scan_id = f"controllerType=0 controllerNumber=1 scan={scan_number}"

    try:
        spectrum = mzml_reader.get_by_id(scan_id)
    except KeyError:
        # Try to find by scan number in ID
        for spec_id in mzml_reader.index.keys():
            if f"scan={scan_number}" in spec_id:
                spectrum = mzml_reader.get_by_id(spec_id)
                break
        else:
            return 'Not found', ''

    # Get filter string from scanList
    filter_string = ''
    if 'scanList' in spectrum and 'scan' in spectrum['scanList']:
        scan_info = spectrum['scanList']['scan'][0]
        filter_string = scan_info.get('filter string', '')

    # Determine activation type from filter string
    filter_lower = filter_string.lower()
    if '@hcd' in filter_lower:
        return 'HCD', filter_string
    elif '@ethcd' in filter_lower:
        return 'EThcD', filter_string
    elif '@etd' in filter_lower:
        return 'ETD', filter_string
    elif '@cid' in filter_lower:
        return 'CID', filter_string
    else:
        # Fallback: check activation element
        if 'precursorList' in spectrum:
            precursor = spectrum['precursorList']['precursor'][0]
            activation = precursor.get('activation', {})
            if 'beam-type collision-induced dissociation' in activation:
                return 'HCD', filter_string
            elif 'electron transfer dissociation' in activation:
                return 'ETD', filter_string

    return 'Unknown', filter_string


# =============================================================================
# Main extraction
# =============================================================================

def main():
    print("=" * 70)
    print("Extracting Activation Types for Low Probability PSMs")
    print("=" * 70)

    # Load low probability PSM file
    if not os.path.exists(LOW_PROB_PSM_FILE):
        print(f"Error: Input file not found: {LOW_PROB_PSM_FILE}")
        return

    df = pd.read_csv(LOW_PROB_PSM_FILE)
    print(f"\nLoaded {len(df)} low probability PSMs")

    # Group PSMs by cell type and mzML file
    file_groups = defaultdict(lambda: defaultdict(list))
    for idx, row in df.iterrows():
        cell_type = row['Cell_Type']
        spectrum_str = row['Spectrum']
        file_name, scan_number, _ = parse_spectrum_id(spectrum_str)
        if file_name and scan_number:
            file_groups[cell_type][file_name].append((idx, scan_number))

    # Results storage
    activation_types = [''] * len(df)
    filter_strings = [''] * len(df)

    # Process each cell type
    for cell_type in ['HEK293T', 'HepG2', 'Jurkat']:
        if cell_type not in file_groups:
            continue

        mzml_dir = MZML_DIRS[cell_type]
        files = file_groups[cell_type]

        print(f"\n{cell_type}: {sum(len(scans) for scans in files.values())} PSMs from {len(files)} files")

        for file_name, scans in files.items():
            mzml_path = find_calibrated_mzml(file_name, mzml_dir)

            if mzml_path is None:
                print(f"  WARNING: Calibrated mzML not found for: {file_name}")
                continue

            print(f"  Processing: {os.path.basename(mzml_path)} ({len(scans)} scans)...", end=" ", flush=True)

            try:
                with mzml.MzML(mzml_path, use_index=True) as reader:
                    extracted = 0
                    for idx, scan_number in scans:
                        act_type, filter_str = get_activation_type(reader, scan_number)
                        activation_types[idx] = act_type
                        filter_strings[idx] = filter_str
                        extracted += 1

                    print(f"done ({extracted})")

            except Exception as e:
                print(f"Error: {e}")

    # Add results to dataframe
    df['Activation_Type'] = activation_types
    df['Filter_String'] = filter_strings

    # Save results
    df.to_csv(OUTPUT_FILE, index=False)
    print(f"\nResults saved to: {OUTPUT_FILE}")

    # =============================================================================
    # Summary
    # =============================================================================

    print("\n" + "=" * 70)
    print("SUMMARY: Activation Type Distribution")
    print("=" * 70)

    # Overall summary
    print("\nOverall:")
    overall_counts = Counter(df['Activation_Type'])
    for act_type, count in overall_counts.most_common():
        pct = count / len(df) * 100
        print(f"  {act_type}: {count} ({pct:.1f}%)")

    # By cell type
    print("\nBy Cell Type:")
    for cell_type in ['HEK293T', 'HepG2', 'Jurkat']:
        cell_data = df[df['Cell_Type'] == cell_type]
        if len(cell_data) > 0:
            print(f"\n  {cell_type} (n={len(cell_data)}):")
            cell_counts = Counter(cell_data['Activation_Type'])
            for act_type, count in cell_counts.most_common():
                pct = count / len(cell_data) * 100
                print(f"    {act_type}: {count} ({pct:.1f}%)")

    # By confidence level
    print("\nBy Confidence Level:")
    for conf in df['Confidence.Level'].unique():
        conf_data = df[df['Confidence.Level'] == conf]
        if len(conf_data) > 0:
            print(f"\n  {conf} (n={len(conf_data)}):")
            conf_counts = Counter(conf_data['Activation_Type'])
            for act_type, count in conf_counts.most_common():
                pct = count / len(conf_data) * 100
                print(f"    {act_type}: {count} ({pct:.1f}%)")

    print("\nDone!")


if __name__ == "__main__":
    main()
