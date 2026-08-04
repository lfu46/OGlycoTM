#!/usr/bin/env python3
"""
Annotate Low Probability Level1b Spectra

This script extracts and annotates all Level1b PSMs with site probability < 0.75
from calibrated mzML files. These are HCD spectra, so we focus on b/y ions and
Y ions (glycopeptide intact ions) for localization evidence analysis.

Output:
- Individual annotated spectrum PDFs
- Combined PDF with all spectra
- Statistics CSV with localization evidence analysis
"""

import os
import sys
import pandas as pd
import numpy as np
from pyteomics import mzml
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import json
import re

# Add current directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fragment_calculator import (
    FragmentCalculator,
    TheoreticalIon,
    match_peaks,
    parse_modifications_from_string,
    calculate_false_match_rate,
    calculate_annotation_statistics
)
from spectrum_annotator import SpectrumAnnotator

# =============================================================================
# Configuration
# =============================================================================

# Input files
LOW_PROB_PSM_FILE = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/site/OGlcNAc_low_prob_psm.csv"

# Site data files (for getting modification info)
SITE_FILES = {
    'HEK293T': '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/site/OGlcNAc_site_HEK293T.csv',
    'HepG2': '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/site/OGlcNAc_site_HepG2.csv',
    'Jurkat': '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/site/OGlcNAc_site_Jurkat.csv'
}

# mzML directories (using calibrated files)
MZML_DIRS = {
    'HEK293T': '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_HEK293T/',
    'HepG2': '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_HepG2/',
    'Jurkat': '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_Jurkat/'
}

# Output directory
OUTPUT_DIR = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/low_prob_spectra_analysis/"


# =============================================================================
# Helper Functions
# =============================================================================

def parse_spectrum_id(spectrum_str):
    """Parse spectrum identifier to extract file name and scan number."""
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


def extract_spectrum_data(mzml_reader, scan_number):
    """Extract spectrum data for a specific scan using indexed access."""
    scan_id = f"controllerType=0 controllerNumber=1 scan={scan_number}"

    try:
        spectrum = mzml_reader.get_by_id(scan_id)
    except KeyError:
        return None

    result = {
        'mz_array': spectrum.get('m/z array', np.array([])),
        'intensity_array': spectrum.get('intensity array', np.array([])),
        'ms_level': spectrum.get('ms level'),
        'tic': spectrum.get('total ion current'),
    }

    # Get precursor info
    if 'precursorList' in spectrum:
        precursor = spectrum['precursorList']['precursor'][0]
        if 'selectedIonList' in precursor:
            sel_ion = precursor['selectedIonList']['selectedIon'][0]
            result['precursor_mz'] = sel_ion.get('selected ion m/z')
            result['precursor_charge'] = sel_ion.get('charge state')

    # Get filter string
    if 'scanList' in spectrum and 'scan' in spectrum['scanList']:
        scan_info = spectrum['scanList']['scan'][0]
        result['filter_string'] = scan_info.get('filter string', '')

    return result


def parse_site_probability(prob_str):
    """Parse site probability string like '[5,I1V1,0.823]'."""
    try:
        match = re.search(r'\[(\d+),([^,]+),([\d.]+)\]', str(prob_str))
        if match:
            return {
                'position': int(match.group(1)),
                'ion_type': match.group(2),
                'probability': float(match.group(3))
            }
    except:
        pass
    return None


def analyze_localization_evidence(annotator, glycan_position, peptide_length):
    """
    Analyze if the spectrum provides evidence for site localization.

    For HCD spectra, localization evidence comes from:
    1. b/y ions that bracket the modification site
    2. Glycan neutral loss ions that indicate modification position

    Returns a dict with localization analysis.
    """
    coverage = annotator._get_fragmentation_coverage()

    # Get all matched ions
    b_ions = coverage.get('b', set())
    y_ions = coverage.get('y', set())

    # Check for bracketing ions
    # b ions: if we have b(n-1) and b(n) where n is the glycan position,
    #         it indicates the modification is at position n
    # y ions: similar logic from C-terminus

    # For N-terminal bracket: need b ions before and after the site
    b_before = any(i < glycan_position for i in b_ions) if glycan_position else False
    b_after = any(i >= glycan_position for i in b_ions) if glycan_position else False
    b_bracket = b_before and b_after

    # For C-terminal bracket: y ions
    y_pos_from_c = peptide_length - glycan_position if glycan_position else 0
    y_before = any(i <= y_pos_from_c for i in y_ions) if glycan_position else False
    y_after = any(i > y_pos_from_c for i in y_ions) if glycan_position else False
    y_bracket = y_before and y_after

    # Count Y ions (glycopeptide intact ions)
    y0_ions = sum(1 for ion in annotator.matched_ions if ion.ion_type == 'Y' and ion.ion_number == 0)
    y1_ions = sum(1 for ion in annotator.matched_ions if ion.ion_type == 'Y' and ion.ion_number == 1)

    # Count oxonium ions (diagnostic for glycan presence)
    oxonium_ions = sum(1 for ion in annotator.matched_ions if ion.ion_type == 'oxonium')

    # Determine localization confidence
    has_bracketing = b_bracket or y_bracket
    has_y_ions = y0_ions > 0 or y1_ions > 0
    has_good_coverage = len(b_ions) + len(y_ions) >= peptide_length * 0.3

    if has_bracketing and has_good_coverage:
        localization_confidence = 'Good'
    elif has_bracketing or has_good_coverage:
        localization_confidence = 'Moderate'
    else:
        localization_confidence = 'Poor'

    return {
        'glycan_position': glycan_position,
        'b_ions_count': len(b_ions),
        'y_ions_count': len(y_ions),
        'b_bracket': b_bracket,
        'y_bracket': y_bracket,
        'has_bracketing': has_bracketing,
        'y0_ions': y0_ions,
        'y1_ions': y1_ions,
        'oxonium_ions': oxonium_ions,
        'localization_confidence': localization_confidence
    }


# =============================================================================
# Main Processing
# =============================================================================

def main():
    print("=" * 70)
    print("Annotating Low Probability Level1b Spectra")
    print("=" * 70)

    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    spectra_dir = os.path.join(OUTPUT_DIR, "individual_spectra")
    os.makedirs(spectra_dir, exist_ok=True)

    # Load low probability PSM file
    if not os.path.exists(LOW_PROB_PSM_FILE):
        print(f"Error: Input file not found: {LOW_PROB_PSM_FILE}")
        return

    low_prob_df = pd.read_csv(LOW_PROB_PSM_FILE)
    print(f"\nLoaded {len(low_prob_df)} low probability PSMs")

    # Load site data files for modification info
    site_data = {}
    for cell_type, file_path in SITE_FILES.items():
        df = pd.read_csv(file_path)
        # Create lookup by spectrum
        for _, row in df.iterrows():
            site_data[row['Spectrum']] = row
    print(f"Loaded site data for {len(site_data)} spectra")

    # Group PSMs by cell type and mzML file
    file_groups = defaultdict(lambda: defaultdict(list))
    for idx, row in low_prob_df.iterrows():
        cell_type = row['Cell_Type']
        spectrum_str = row['Spectrum']
        file_name, scan_number, _ = parse_spectrum_id(spectrum_str)
        if file_name and scan_number:
            file_groups[cell_type][file_name].append({
                'idx': idx,
                'scan': scan_number,
                'row': row,
                'spectrum': spectrum_str
            })

    # Results storage
    results = []
    annotated_count = 0
    failed_count = 0

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
                failed_count += len(scans)
                continue

            print(f"  Processing: {os.path.basename(mzml_path)} ({len(scans)} scans)...", end=" ", flush=True)

            try:
                with mzml.MzML(mzml_path, use_index=True) as reader:
                    extracted = 0

                    for scan_info in scans:
                        idx = scan_info['idx']
                        scan_number = scan_info['scan']
                        row = scan_info['row']
                        spectrum_str = scan_info['spectrum']

                        # Extract spectrum
                        spec_data = extract_spectrum_data(reader, scan_number)
                        if spec_data is None:
                            failed_count += 1
                            continue

                        # Get full site data for modifications
                        if spectrum_str in site_data:
                            full_row = site_data[spectrum_str]
                            assigned_mods = full_row.get('Assigned.Modifications', '')
                            peptide = full_row['Peptide']
                            charge = int(full_row['Charge'])
                            precursor_mz = float(full_row.get('Observed.M.Z', spec_data.get('precursor_mz', 0)))
                        else:
                            assigned_mods = ''
                            peptide = row['Peptide']
                            charge = int(row['Charge'])
                            precursor_mz = spec_data.get('precursor_mz', 0)

                        # Parse modifications
                        modifications = parse_modifications_from_string(assigned_mods)

                        # Parse site probability info
                        prob_info = parse_site_probability(row.get('Site.Probabilities', ''))
                        site_prob = row.get('Site_Prob', prob_info['probability'] if prob_info else 0)

                        try:
                            # Create annotator
                            annotator = SpectrumAnnotator(
                                peptide=peptide,
                                modifications=modifications,
                                precursor_charge=charge,
                                precursor_mz=precursor_mz,
                                exp_mz=spec_data['mz_array'],
                                exp_intensity=spec_data['intensity_array'],
                                tolerance_ppm=20.0,
                                site_index=row['site_index'],
                                gene=row['Gene']
                            )

                            # Get glycan position
                            glycan_position = annotator.glycan_position

                            # Analyze localization evidence
                            loc_evidence = analyze_localization_evidence(
                                annotator, glycan_position, len(peptide)
                            )

                            # Generate output filename
                            output_file = os.path.join(
                                spectra_dir,
                                f"{cell_type}_{row['site_index']}_{scan_number}.pdf"
                            )

                            # Create plot
                            fig = annotator.plot(output_path=output_file, show_error_plot=True)
                            plt.close(fig)

                            # Collect statistics
                            fmr = annotator.false_match_rate
                            stats = annotator.annotation_stats

                            results.append({
                                'cell_type': cell_type,
                                'site_index': row['site_index'],
                                'gene': row['Gene'],
                                'peptide': peptide,
                                'charge': charge,
                                'scan_number': scan_number,
                                'confidence_level': row['Confidence.Level'],
                                'site_probability': site_prob,
                                'glycan_position': glycan_position,
                                'sequence_coverage': stats['sequence_coverage'],
                                'intensity_annotated': stats['intensity_annotated'],
                                'fmr_peaks': fmr.fmr_peaks,
                                'fmr_intensity': fmr.fmr_intensity,
                                'b_ions_count': loc_evidence['b_ions_count'],
                                'y_ions_count': loc_evidence['y_ions_count'],
                                'b_bracket': loc_evidence['b_bracket'],
                                'y_bracket': loc_evidence['y_bracket'],
                                'has_bracketing': loc_evidence['has_bracketing'],
                                'y0_ions': loc_evidence['y0_ions'],
                                'y1_ions': loc_evidence['y1_ions'],
                                'oxonium_ions': loc_evidence['oxonium_ions'],
                                'localization_confidence': loc_evidence['localization_confidence'],
                                'output_file': output_file
                            })

                            annotated_count += 1
                            extracted += 1

                        except Exception as e:
                            print(f"\n    Error annotating {row['site_index']}: {e}")
                            failed_count += 1
                            continue

                    print(f"done ({extracted})")

            except Exception as e:
                print(f"Error: {e}")
                failed_count += len(scans)

    # Save results
    results_df = pd.DataFrame(results)
    stats_file = os.path.join(OUTPUT_DIR, "low_prob_spectra_analysis.csv")
    results_df.to_csv(stats_file, index=False)
    print(f"\nSaved analysis to: {stats_file}")

    # ==========================================================================
    # Combine PDFs
    # ==========================================================================
    print("\nCombining PDFs...")
    combined_pdf_path = os.path.join(OUTPUT_DIR, "all_low_prob_spectra_annotated.pdf")

    # Get all PDF files sorted by cell type and site
    pdf_files = sorted(results_df['output_file'].tolist())

    if pdf_files:
        from PyPDF2 import PdfMerger

        merger = PdfMerger()
        for pdf_file in pdf_files:
            if os.path.exists(pdf_file):
                merger.append(pdf_file)

        merger.write(combined_pdf_path)
        merger.close()
        print(f"Combined PDF saved to: {combined_pdf_path}")
    else:
        print("No PDFs to combine")

    # ==========================================================================
    # Summary Statistics
    # ==========================================================================
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    print(f"\nTotal PSMs processed: {annotated_count}")
    print(f"Failed: {failed_count}")

    if len(results_df) > 0:
        print("\nLocalization Evidence Analysis:")
        print("-" * 50)

        # By localization confidence
        conf_counts = results_df['localization_confidence'].value_counts()
        print("\nLocalization Confidence:")
        for conf, count in conf_counts.items():
            pct = count / len(results_df) * 100
            print(f"  {conf}: {count} ({pct:.1f}%)")

        # By cell type
        print("\nBy Cell Type:")
        for cell_type in ['HEK293T', 'HepG2', 'Jurkat']:
            cell_data = results_df[results_df['cell_type'] == cell_type]
            if len(cell_data) > 0:
                good = sum(cell_data['localization_confidence'] == 'Good')
                moderate = sum(cell_data['localization_confidence'] == 'Moderate')
                poor = sum(cell_data['localization_confidence'] == 'Poor')
                print(f"  {cell_type}: Good={good}, Moderate={moderate}, Poor={poor}")

        # Average statistics
        print("\nAverage Statistics:")
        print(f"  Sequence Coverage: {results_df['sequence_coverage'].mean()*100:.1f}%")
        print(f"  Intensity Annotated: {results_df['intensity_annotated'].mean()*100:.1f}%")
        print(f"  FMR (peaks): {results_df['fmr_peaks'].mean()*100:.1f}%")
        print(f"  Has Bracketing Ions: {results_df['has_bracketing'].mean()*100:.1f}%")

    print("\n" + "=" * 70)
    print("Done!")
    print("=" * 70)


if __name__ == "__main__":
    main()
