#!/usr/bin/env python3
"""
Analyze Y1 Ion Detection Across All EThcD Spectra

Y1 ions represent the intact glycopeptide (peptide + full glycan).
This script analyzes how often Y1 ions are detected and their intensities.
"""

import os
import sys
import json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from fragment_calculator import (
    FragmentCalculator,
    match_peaks,
    PROTON,
)

SOURCE_PATH = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/"

def analyze_Y1_ions(cell_type: str, tolerance_ppm: float = 20.0):
    """Analyze Y1 ion detection for a cell type."""

    summary_file = os.path.join(
        SOURCE_PATH,
        f"point_to_point_response/extracted_spectra_EThcD_{cell_type}_summary.csv"
    )
    spectra_dir = os.path.join(
        SOURCE_PATH,
        f"point_to_point_response/extracted_spectra_EThcD_{cell_type}"
    )

    df = pd.read_csv(summary_file)

    results = []

    for idx, row in df.iterrows():
        try:
            # Load spectrum
            spec_file = os.path.join(spectra_dir, row['spectrum_file'])
            spec_df = pd.read_csv(spec_file)
            exp_mz = spec_df['mz'].values
            exp_intensity = spec_df['intensity'].values

            # Parse modifications
            if 'modifications_json' in row and pd.notna(row['modifications_json']):
                modifications = json.loads(row['modifications_json'])
            else:
                modifications = []

            # Create calculator
            calculator = FragmentCalculator(
                row['Peptide'],
                modifications,
                int(row['Charge']),
                max_fragment_charge=2
            )

            # Get Y ions
            y_ions = calculator.calculate_Y_ions()

            # Separate Y0 and Y1 ions
            y0_ions = [ion for ion in y_ions if ion.ion_number == 0]
            y1_ions = [ion for ion in y_ions if ion.ion_number == 1]

            # Match Y1 ions
            base_peak = np.max(exp_intensity)

            y1_matched = []
            for ion in y1_ions:
                # Find matching peak
                ppm_errors = np.abs((exp_mz - ion.mz) / ion.mz * 1e6)
                matches = np.where(ppm_errors <= tolerance_ppm)[0]

                if len(matches) > 0:
                    best_match = matches[np.argmax(exp_intensity[matches])]
                    rel_intensity = exp_intensity[best_match] / base_peak * 100
                    y1_matched.append({
                        'charge': ion.charge,
                        'theoretical_mz': ion.mz,
                        'exp_mz': exp_mz[best_match],
                        'intensity': exp_intensity[best_match],
                        'rel_intensity': rel_intensity,
                        'ppm_error': ppm_errors[best_match]
                    })

            # Match Y0 ions
            y0_matched = []
            for ion in y0_ions:
                ppm_errors = np.abs((exp_mz - ion.mz) / ion.mz * 1e6)
                matches = np.where(ppm_errors <= tolerance_ppm)[0]

                if len(matches) > 0:
                    best_match = matches[np.argmax(exp_intensity[matches])]
                    rel_intensity = exp_intensity[best_match] / base_peak * 100
                    y0_matched.append({
                        'charge': ion.charge,
                        'theoretical_mz': ion.mz,
                        'exp_mz': exp_mz[best_match],
                        'intensity': exp_intensity[best_match],
                        'rel_intensity': rel_intensity,
                        'ppm_error': ppm_errors[best_match]
                    })

            results.append({
                'site_index': row['site_index'],
                'scan_number': row['scan_number'],
                'gene': row['Gene'],
                'peptide': row['Peptide'],
                'charge': int(row['Charge']),
                'n_y1_theoretical': len(y1_ions),
                'n_y1_matched': len(y1_matched),
                'y1_detected': len(y1_matched) > 0,
                'y1_max_rel_intensity': max([m['rel_intensity'] for m in y1_matched]) if y1_matched else 0,
                'y1_charges_detected': [m['charge'] for m in y1_matched],
                'n_y0_theoretical': len(y0_ions),
                'n_y0_matched': len(y0_matched),
                'y0_detected': len(y0_matched) > 0,
                'y0_max_rel_intensity': max([m['rel_intensity'] for m in y0_matched]) if y0_matched else 0,
            })

        except Exception as e:
            print(f"  Error processing {row['site_index']}: {e}")
            continue

    return pd.DataFrame(results)


def main():
    print("=" * 70)
    print("Y1 Ion Analysis Across All EThcD Spectra")
    print("=" * 70)
    print("\nY1 = Intact glycopeptide (peptide + full glycan)")
    print("Y0 = Peptide only (glycan lost)")
    print()

    all_results = []

    for cell_type in ['HEK293T', 'HepG2', 'Jurkat']:
        print(f"\n{'='*50}")
        print(f"Analyzing {cell_type}...")
        print('='*50)

        results = analyze_Y1_ions(cell_type)
        results['cell_type'] = cell_type
        all_results.append(results)

        # Summary statistics
        n_total = len(results)
        n_y1_detected = results['y1_detected'].sum()
        n_y0_detected = results['y0_detected'].sum()

        print(f"\nTotal spectra: {n_total}")
        print(f"\nY1 (intact glycopeptide) detection:")
        print(f"  Spectra with Y1 detected: {n_y1_detected}/{n_total} ({n_y1_detected/n_total*100:.1f}%)")

        y1_intensities = results[results['y1_detected']]['y1_max_rel_intensity']
        if len(y1_intensities) > 0:
            print(f"  Y1 relative intensity (% base peak):")
            print(f"    Mean: {y1_intensities.mean():.1f}%")
            print(f"    Median: {y1_intensities.median():.1f}%")
            print(f"    Range: {y1_intensities.min():.1f}% - {y1_intensities.max():.1f}%")

        print(f"\nY0 (glycan loss) detection:")
        print(f"  Spectra with Y0 detected: {n_y0_detected}/{n_total} ({n_y0_detected/n_total*100:.1f}%)")

        y0_intensities = results[results['y0_detected']]['y0_max_rel_intensity']
        if len(y0_intensities) > 0:
            print(f"  Y0 relative intensity (% base peak):")
            print(f"    Mean: {y0_intensities.mean():.1f}%")
            print(f"    Median: {y0_intensities.median():.1f}%")
            print(f"    Range: {y0_intensities.min():.1f}% - {y0_intensities.max():.1f}%")

    # Combined results
    combined = pd.concat(all_results, ignore_index=True)

    print("\n" + "=" * 70)
    print("OVERALL SUMMARY (All Cell Types)")
    print("=" * 70)

    n_total = len(combined)
    n_y1_detected = combined['y1_detected'].sum()
    n_y0_detected = combined['y0_detected'].sum()

    print(f"\nTotal spectra: {n_total}")
    print(f"\nY1 (intact glycopeptide):")
    print(f"  Detected in: {n_y1_detected}/{n_total} spectra ({n_y1_detected/n_total*100:.1f}%)")

    y1_intensities = combined[combined['y1_detected']]['y1_max_rel_intensity']
    if len(y1_intensities) > 0:
        print(f"  Relative intensity when detected:")
        print(f"    Mean: {y1_intensities.mean():.1f}%")
        print(f"    Median: {y1_intensities.median():.1f}%")
        print(f"    25th percentile: {y1_intensities.quantile(0.25):.1f}%")
        print(f"    75th percentile: {y1_intensities.quantile(0.75):.1f}%")

    print(f"\nY0 (glycan loss):")
    print(f"  Detected in: {n_y0_detected}/{n_total} spectra ({n_y0_detected/n_total*100:.1f}%)")

    y0_intensities = combined[combined['y0_detected']]['y0_max_rel_intensity']
    if len(y0_intensities) > 0:
        print(f"  Relative intensity when detected:")
        print(f"    Mean: {y0_intensities.mean():.1f}%")
        print(f"    Median: {y0_intensities.median():.1f}%")

    # Save detailed results
    output_file = os.path.join(SOURCE_PATH, "point_to_point_response/Y1_ion_analysis.csv")
    combined.to_csv(output_file, index=False)
    print(f"\nDetailed results saved to: {output_file}")

    print("\n" + "=" * 70)

    return combined


if __name__ == "__main__":
    results = main()
