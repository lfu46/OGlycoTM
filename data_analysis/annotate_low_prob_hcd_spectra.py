#!/usr/bin/env python3
"""
Annotate Low Probability HCD Spectra for O-GlcNAc Site Localization Analysis

This script uses the GlycoSpectrumAnnotator package to annotate all 471 Level1b PSMs
with site probability < 0.75. These are HCD spectra where site localization
confidence is limited compared to EThcD spectra.

The goal is to:
1. Extract spectra from calibrated mzML files
2. Annotate each spectrum with theoretical fragment ions
3. Calculate annotation statistics and False Match Rate
4. Analyze localization evidence (site-determining ions)
5. Combine results into a PDF report

Author: Claude Code / Longping Fu
"""

import os
import sys
import json
import numpy as np
import pandas as pd
from pyteomics import mzml
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Add GlycoSpectrumAnnotator to path
sys.path.insert(0, '/Users/longpingfu/Downloads/GlycoSpectrumAnnotator')

from spectrum_annotator_ddzby import (
    SpectrumAnnotator,
    FragmentCalculator,
    parse_modifications_from_string,
    match_peaks,
    calculate_false_match_rate,
    calculate_annotation_statistics,
)

# =============================================================================
# Configuration
# =============================================================================

# Data paths
SOURCE_PATH = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/"
OUTPUT_PATH = os.path.join(SOURCE_PATH, "low_prob_spectra/")

# Input file with activation types
LOW_PROB_FILE = os.path.join(SOURCE_PATH, "site/OGlcNAc_low_prob_psm_activation.csv")

# Original bonafide files (to get Assigned.Modifications)
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


# =============================================================================
# Helper Functions
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


def extract_spectrum_data(mzml_reader, scan_number):
    """
    Extract complete spectrum data for a specific scan using indexed access.
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
            return None

    result = {
        'mz_array': spectrum.get('m/z array', np.array([])),
        'intensity_array': spectrum.get('intensity array', np.array([])),
        'ms_level': spectrum.get('ms level'),
    }

    # Get precursor info
    if 'precursorList' in spectrum:
        precursor = spectrum['precursorList']['precursor'][0]
        if 'selectedIonList' in precursor:
            sel_ion = precursor['selectedIonList']['selectedIon'][0]
            result['precursor_mz'] = sel_ion.get('selected ion m/z')
            result['precursor_charge'] = sel_ion.get('charge state')

    return result


def analyze_site_determining_ions(peptide, glycan_position, matched_ions, peptide_length):
    """
    Analyze which matched ions provide evidence for site localization.

    Site-determining ions are fragment ions that:
    - Contain only one potential glycosylation site (S/T/Y)
    - Have or lack the glycan mass depending on which site is modified

    For HCD spectra, these are primarily b and y ions.
    """
    # Find all potential sites (S, T, Y positions)
    potential_sites = []
    for i, aa in enumerate(peptide):
        if aa in ['S', 'T', 'Y']:
            potential_sites.append(i + 1)  # 1-indexed

    if len(potential_sites) <= 1:
        # Only one possible site - automatically localized
        return {
            'n_potential_sites': len(potential_sites),
            'site_determining_ions': [],
            'localization_evidence': 'single_site'
        }

    # Analyze matched b and y ions for site localization
    site_determining_ions = []

    for ion in matched_ions:
        if ion.ion_type not in ['b', 'y']:
            continue

        # Determine which sites are covered by this ion
        if ion.ion_type == 'b':
            # b_n contains residues 1 to n
            covered_sites = [s for s in potential_sites if s <= ion.ion_number]
        else:  # y ion
            # y_n contains residues (L-n+1) to L
            first_residue = peptide_length - ion.ion_number + 1
            covered_sites = [s for s in potential_sites if s >= first_residue]

        # If ion covers exactly one site, it's site-determining
        if len(covered_sites) == 1:
            site_determining_ions.append({
                'ion': ion.annotation,
                'exp_mz': ion.exp_mz,
                'theo_mz': ion.mz,
                'covered_site': covered_sites[0],
                'has_glycan': glycan_position in covered_sites,
                'intensity': ion.exp_intensity
            })

    return {
        'n_potential_sites': len(potential_sites),
        'potential_sites': potential_sites,
        'glycan_position': glycan_position,
        'site_determining_ions': site_determining_ions,
        'n_site_determining_ions': len(site_determining_ions),
        'localization_evidence': 'multiple_sites'
    }


# =============================================================================
# Main Processing
# =============================================================================

def main():
    print("=" * 70)
    print("Annotating Low Probability HCD Spectra")
    print("Using GlycoSpectrumAnnotator Package")
    print("=" * 70)

    # Create output directory
    os.makedirs(OUTPUT_PATH, exist_ok=True)
    spectra_dir = os.path.join(OUTPUT_PATH, "annotated_spectra")
    os.makedirs(spectra_dir, exist_ok=True)

    # Load low probability PSM file
    if not os.path.exists(LOW_PROB_FILE):
        print(f"Error: Input file not found: {LOW_PROB_FILE}")
        return

    df_low_prob = pd.read_csv(LOW_PROB_FILE)
    print(f"\nLoaded {len(df_low_prob)} low probability PSMs")

    # Load bonafide files to get Assigned.Modifications
    bonafide_dfs = {}
    for cell_type, path in BONAFIDE_FILES.items():
        if os.path.exists(path):
            bonafide_dfs[cell_type] = pd.read_csv(path)
            print(f"  Loaded {len(bonafide_dfs[cell_type])} bonafide PSMs for {cell_type}")

    # Merge to get Assigned.Modifications
    merged_data = []
    for idx, row in df_low_prob.iterrows():
        cell_type = row['Cell_Type']
        spectrum = row['Spectrum']

        if cell_type in bonafide_dfs:
            # Find matching PSM
            match = bonafide_dfs[cell_type][bonafide_dfs[cell_type]['Spectrum'] == spectrum]
            if len(match) > 0:
                merged_row = row.to_dict()
                merged_row['Assigned.Modifications'] = match.iloc[0].get('Assigned.Modifications', '')
                merged_data.append(merged_row)

    df_merged = pd.DataFrame(merged_data)
    print(f"\nMerged data: {len(df_merged)} PSMs with modification info")

    # Group by cell type and mzML file for efficient processing
    file_groups = defaultdict(lambda: defaultdict(list))
    for idx, row in df_merged.iterrows():
        cell_type = row['Cell_Type']
        spectrum_str = row['Spectrum']
        file_name, scan_number, _ = parse_spectrum_id(spectrum_str)
        if file_name and scan_number:
            file_groups[cell_type][file_name].append((idx, scan_number, row))

    # Results storage
    results = []
    n_processed = 0
    n_errors = 0

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

                        # Parse modifications
                        mod_string = row.get('Assigned.Modifications', '')
                        modifications = parse_modifications_from_string(mod_string if pd.notna(mod_string) else '')

                        # Get precursor info
                        precursor_mz = row.get('Observed.M.Z', spec_data.get('precursor_mz', 0))
                        if pd.isna(precursor_mz):
                            precursor_mz = spec_data.get('precursor_mz', 0)
                        precursor_charge = int(row['Charge'])

                        peptide = row['Peptide']

                        # Find glycan position from modifications
                        glycan_position = None
                        for mod in modifications:
                            if abs(mod['mass'] - 528.2859) < 0.1:  # HexNAc+TMT
                                glycan_position = mod['position']
                                break

                        try:
                            # Create annotator (HCD mode - only b/y ions, no c/z)
                            annotator = SpectrumAnnotator(
                                peptide=peptide,
                                modifications=modifications,
                                precursor_charge=precursor_charge,
                                precursor_mz=float(precursor_mz) if precursor_mz else 0,
                                exp_mz=exp_mz,
                                exp_intensity=exp_intensity,
                                tolerance_ppm=20.0,
                                site_index=row['site_index'],
                                gene=row['Gene'],
                                activation_type="HCD"  # Only b/y/Y/oxonium ions, no c/z
                            )

                            # Analyze site-determining ions
                            site_analysis = analyze_site_determining_ions(
                                peptide, glycan_position, annotator.matched_ions, len(peptide)
                            )

                            # Collect results
                            fmr = annotator.false_match_rate
                            stats = annotator.annotation_stats

                            result = {
                                # Identifiers
                                'Cell_Type': cell_type,
                                'Spectrum': row['Spectrum'],
                                'site_index': row['site_index'],
                                'Gene': row['Gene'],
                                'scan_number': scan_number,

                                # Peptide info
                                'Peptide': peptide,
                                'Charge': precursor_charge,
                                'Assigned.Modifications': mod_string,
                                'Site_Prob': row['Site_Prob'],
                                'Confidence.Level': row['Confidence.Level'],

                                # Annotation statistics
                                'sequence_coverage': stats['sequence_coverage'],
                                'sequence_coverage_bonds': stats['sequence_coverage_bonds'],
                                'peaks_annotated': stats['peaks_annotated'],
                                'peaks_annotated_count': stats['peaks_annotated_count'],
                                'intensity_annotated': stats['intensity_annotated'],

                                # False match rate
                                'fmr_peaks': fmr.fmr_peaks,
                                'fmr_intensity': fmr.fmr_intensity,
                                'matched_peaks': fmr.matched_peaks,

                                # Site localization analysis
                                'n_potential_sites': site_analysis['n_potential_sites'],
                                'n_site_determining_ions': site_analysis.get('n_site_determining_ions', 0),
                                'glycan_position': glycan_position,
                                'site_determining_ions_json': json.dumps(site_analysis.get('site_determining_ions', [])),

                                # Spectrum quality
                                'n_peaks': len(exp_mz),
                                'tic': np.sum(exp_intensity),
                                'base_peak_intensity': np.max(exp_intensity),
                            }

                            results.append(result)

                            # Generate PDF for this spectrum
                            output_file = os.path.join(spectra_dir, f"{row['site_index']}_{scan_number}.pdf")
                            fig = annotator.plot(output_path=output_file)
                            plt.close(fig)

                            scans_processed += 1
                            n_processed += 1

                        except Exception as e:
                            n_errors += 1
                            print(f"\n    Error annotating {row['site_index']}: {e}")
                            continue

                    print(f"done ({scans_processed})")

            except Exception as e:
                print(f"Error opening mzML: {e}")
                n_errors += len(scans)

    # Save results
    results_df = pd.DataFrame(results)
    results_file = os.path.join(OUTPUT_PATH, "low_prob_annotation_results.csv")
    results_df.to_csv(results_file, index=False)
    print(f"\nResults saved to: {results_file}")

    # ==========================================================================
    # Summary Statistics
    # ==========================================================================
    print("\n" + "=" * 70)
    print("ANNOTATION SUMMARY")
    print("=" * 70)

    print(f"\nProcessing:")
    print(f"  Total PSMs: {len(df_low_prob)}")
    print(f"  Successfully annotated: {n_processed}")
    print(f"  Errors: {n_errors}")

    if len(results_df) > 0:
        print(f"\nAnnotation Statistics (mean ± std):")
        print(f"  Sequence Coverage: {results_df['sequence_coverage'].mean()*100:.1f} ± {results_df['sequence_coverage'].std()*100:.1f}%")
        print(f"  Intensity Annotated: {results_df['intensity_annotated'].mean()*100:.1f} ± {results_df['intensity_annotated'].std()*100:.1f}%")
        print(f"  FMR (peaks): {results_df['fmr_peaks'].mean()*100:.1f} ± {results_df['fmr_peaks'].std()*100:.1f}%")

        print(f"\nSite Localization Analysis:")
        print(f"  PSMs with single potential site: {(results_df['n_potential_sites'] == 1).sum()}")
        print(f"  PSMs with multiple potential sites: {(results_df['n_potential_sites'] > 1).sum()}")

        multi_site = results_df[results_df['n_potential_sites'] > 1]
        if len(multi_site) > 0:
            print(f"  Mean site-determining ions (multi-site): {multi_site['n_site_determining_ions'].mean():.1f}")

        # By cell type
        print(f"\nBy Cell Type:")
        for cell_type in ['HEK293T', 'HepG2', 'Jurkat']:
            cell_data = results_df[results_df['Cell_Type'] == cell_type]
            if len(cell_data) > 0:
                print(f"  {cell_type}: {len(cell_data)} spectra, "
                      f"seq cov: {cell_data['sequence_coverage'].mean()*100:.1f}%, "
                      f"FMR: {cell_data['fmr_peaks'].mean()*100:.1f}%")

    # ==========================================================================
    # Generate Combined PDF Report
    # ==========================================================================
    print(f"\nGenerating combined PDF report...")

    combined_pdf_path = os.path.join(OUTPUT_PATH, "low_prob_spectra_annotated.pdf")

    with PdfPages(combined_pdf_path) as pdf:
        # Title page
        fig = plt.figure(figsize=(8.5, 11))
        fig.text(0.5, 0.7, "Low Probability HCD Spectra", ha='center', fontsize=20, fontweight='bold')
        fig.text(0.5, 0.6, "Annotation Analysis", ha='center', fontsize=16)
        fig.text(0.5, 0.5, f"Total Spectra: {n_processed}", ha='center', fontsize=12)
        fig.text(0.5, 0.45, f"Level1b PSMs with Site Probability < 0.75", ha='center', fontsize=10)
        plt.axis('off')
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

        # Summary statistics page
        if len(results_df) > 0:
            fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))

            # Distribution of sequence coverage
            axes[0, 0].hist(results_df['sequence_coverage'] * 100, bins=20, edgecolor='black', alpha=0.7)
            axes[0, 0].set_xlabel('Sequence Coverage (%)')
            axes[0, 0].set_ylabel('Count')
            axes[0, 0].set_title('Sequence Coverage Distribution')
            axes[0, 0].axvline(results_df['sequence_coverage'].mean() * 100, color='red', linestyle='--',
                              label=f'Mean: {results_df["sequence_coverage"].mean()*100:.1f}%')
            axes[0, 0].legend()

            # Distribution of FMR
            axes[0, 1].hist(results_df['fmr_peaks'] * 100, bins=20, edgecolor='black', alpha=0.7, color='orange')
            axes[0, 1].set_xlabel('False Match Rate (%)')
            axes[0, 1].set_ylabel('Count')
            axes[0, 1].set_title('False Match Rate Distribution')
            axes[0, 1].axvline(results_df['fmr_peaks'].mean() * 100, color='red', linestyle='--',
                              label=f'Mean: {results_df["fmr_peaks"].mean()*100:.1f}%')
            axes[0, 1].legend()

            # Site probability vs sequence coverage
            axes[1, 0].scatter(results_df['Site_Prob'], results_df['sequence_coverage'] * 100, alpha=0.5)
            axes[1, 0].set_xlabel('Site Probability')
            axes[1, 0].set_ylabel('Sequence Coverage (%)')
            axes[1, 0].set_title('Site Probability vs Sequence Coverage')

            # Number of site-determining ions
            multi_site = results_df[results_df['n_potential_sites'] > 1]
            if len(multi_site) > 0:
                axes[1, 1].hist(multi_site['n_site_determining_ions'], bins=range(0, 15),
                               edgecolor='black', alpha=0.7, color='green')
                axes[1, 1].set_xlabel('Number of Site-Determining Ions')
                axes[1, 1].set_ylabel('Count')
                axes[1, 1].set_title('Site-Determining Ions (Multi-Site PSMs)')
            else:
                axes[1, 1].text(0.5, 0.5, 'No multi-site PSMs', ha='center', va='center')
                axes[1, 1].axis('off')

            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)

        # Add individual spectrum PDFs
        spectrum_files = sorted([f for f in os.listdir(spectra_dir) if f.endswith('.pdf')])
        print(f"  Adding {len(spectrum_files)} individual spectra...")

        # Note: We can't easily merge PDFs with PdfPages, so we'll note the location
        print(f"\n  Individual spectrum PDFs saved in: {spectra_dir}/")

    print(f"\nCombined PDF saved to: {combined_pdf_path}")
    print("\nDone!")


if __name__ == "__main__":
    main()
