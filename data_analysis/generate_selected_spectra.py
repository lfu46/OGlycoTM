#!/usr/bin/env python3
"""
Generate Selected Spectra in Multiple Formats (PDF, TIFF, EMF)

This script regenerates the 6 selected spectra for the reviewer response
and saves them in multiple high-quality formats.

Output formats:
- PDF: Vector format (original)
- TIFF: 600 DPI raster (publication quality)
- EMF: Enhanced Metafile via SVG conversion (Windows vector format)

Author: Claude Code Assistant
Date: 2026-02-03
"""

import os
import sys
import json
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Add GlycoSpectrumAnnotator to path for updated annotator
sys.path.insert(0, '/Users/longpingfu/Downloads/GlycoSpectrumAnnotator')

from spectrum_annotator_ddzby import (
    SpectrumAnnotator,
    FragmentCalculator,
    parse_modifications_from_string,
)

# =============================================================================
# CONFIGURATION
# =============================================================================

SOURCE_PATH = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/"
OUTPUT_DIR = os.path.join(SOURCE_PATH, "point_to_point_response/selected_spectra")

# Selected spectra to regenerate (site_index, scan_number, cell_type)
SELECTED_SPECTRA = [
    ("O00268_S528", 11775, "Jurkat"),
    ("P49790_S893", 34903, "HepG2"),
    ("P51610_T490", 37826, "HepG2"),
    ("Q14157_S445", 35762, "Jurkat"),
    ("Q8WXI9_S584", 22196, "HepG2"),
    ("Q9NYV4_S593", 32543, "HepG2"),
]

# Output settings
TIFF_DPI = 600
FIGSIZE = (7, 5)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def load_spectrum_metadata(site_index: str, scan_number: int, cell_type: str) -> dict:
    """Load metadata for a spectrum from the summary CSV."""
    summary_file = os.path.join(
        SOURCE_PATH,
        f"point_to_point_response/extracted_spectra_EThcD_{cell_type}_summary.csv"
    )

    df = pd.read_csv(summary_file)

    # Find the matching row
    mask = (df['site_index'] == site_index) & (df['scan_number'] == scan_number)
    if mask.sum() == 0:
        raise ValueError(f"Spectrum not found: {site_index}_{scan_number} in {cell_type}")

    return df[mask].iloc[0].to_dict()


def load_spectrum_data(site_index: str, scan_number: int, cell_type: str) -> tuple:
    """Load m/z and intensity arrays for a spectrum."""
    spectrum_file = os.path.join(
        SOURCE_PATH,
        f"point_to_point_response/extracted_spectra_EThcD_{cell_type}",
        f"{site_index}_{scan_number}.csv"
    )

    df = pd.read_csv(spectrum_file)
    return df['mz'].values, df['intensity'].values


def svg_to_emf(svg_path: str, emf_path: str) -> bool:
    """Convert SVG to EMF using Inkscape."""
    inkscape_path = "/opt/homebrew/bin/inkscape"

    if not os.path.exists(inkscape_path):
        print(f"  Warning: Inkscape not found at {inkscape_path}")
        return False

    try:
        result = subprocess.run(
            [inkscape_path, svg_path, f"--export-filename={emf_path}"],
            capture_output=True,
            text=True,
            timeout=60
        )
        if result.returncode == 0:
            return True
        else:
            print(f"  Warning: Inkscape conversion failed: {result.stderr[:200]}")
            return False
    except Exception as e:
        print(f"  Warning: SVG to EMF conversion error: {e}")
        return False


# =============================================================================
# MAIN FUNCTION
# =============================================================================

def generate_spectrum(site_index: str, scan_number: int, cell_type: str,
                     output_dir: str, tolerance_ppm: float = 20.0):
    """
    Generate annotated spectrum and save in multiple formats.

    Args:
        site_index: Site identifier (e.g., "O00268_S528")
        scan_number: MS scan number
        cell_type: Cell type (HEK293T, HepG2, Jurkat)
        output_dir: Output directory
        tolerance_ppm: Mass tolerance for peak matching

    Returns:
        dict with output file paths
    """
    print(f"\nProcessing: {site_index}_{scan_number} ({cell_type})")

    # Load metadata
    meta = load_spectrum_metadata(site_index, scan_number, cell_type)

    # Load spectrum data
    exp_mz, exp_intensity = load_spectrum_data(site_index, scan_number, cell_type)

    # Parse modifications
    if 'modifications_json' in meta and pd.notna(meta['modifications_json']):
        modifications = json.loads(meta['modifications_json'])
    else:
        modifications = []

    # Create annotator (EThcD spectra - use all ion types including c/z)
    annotator = SpectrumAnnotator(
        peptide=meta['Peptide'],
        modifications=modifications,
        precursor_charge=int(meta['Charge']),
        precursor_mz=float(meta.get('Calibrated_Observed_MZ', meta.get('mzML_precursor_mz', 0))),
        exp_mz=exp_mz,
        exp_intensity=exp_intensity,
        tolerance_ppm=tolerance_ppm,
        site_index=meta['site_index'],
        gene=meta['Gene'],
        activation_type="EThcD"
    )

    # Base filename
    base_name = f"{site_index}_{scan_number}"

    # Create output paths
    pdf_path = os.path.join(output_dir, f"{base_name}.pdf")
    tiff_path = os.path.join(output_dir, f"{base_name}.tiff")
    svg_path = os.path.join(output_dir, f"{base_name}.svg")
    emf_path = os.path.join(output_dir, f"{base_name}.emf")

    # Generate figure
    fig = annotator.plot(figsize=FIGSIZE, output_path=None)

    # Save PDF
    fig.savefig(pdf_path, format='pdf', dpi=300, bbox_inches='tight')
    print(f"  Saved PDF: {pdf_path}")

    # Save TIFF (600 DPI)
    fig.savefig(tiff_path, format='tiff', dpi=TIFF_DPI, bbox_inches='tight')
    print(f"  Saved TIFF: {tiff_path} ({TIFF_DPI} DPI)")

    # Save SVG (intermediate for EMF conversion)
    fig.savefig(svg_path, format='svg', bbox_inches='tight')
    print(f"  Saved SVG: {svg_path}")

    # Convert SVG to EMF
    if svg_to_emf(svg_path, emf_path):
        print(f"  Saved EMF: {emf_path}")
        # Optionally remove SVG after successful EMF conversion
        # os.remove(svg_path)
    else:
        print(f"  EMF conversion failed, keeping SVG")

    plt.close(fig)

    return {
        'pdf': pdf_path,
        'tiff': tiff_path,
        'svg': svg_path,
        'emf': emf_path if os.path.exists(emf_path) else None
    }


def main():
    """Main function to regenerate all selected spectra."""
    print("=" * 70)
    print("Generating Selected Spectra in Multiple Formats")
    print("=" * 70)
    print(f"\nOutput directory: {OUTPUT_DIR}")
    print(f"TIFF resolution: {TIFF_DPI} DPI")
    print(f"Number of spectra: {len(SELECTED_SPECTRA)}")

    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Process each spectrum
    results = []
    for site_index, scan_number, cell_type in SELECTED_SPECTRA:
        try:
            output_files = generate_spectrum(
                site_index, scan_number, cell_type, OUTPUT_DIR
            )
            results.append({
                'site_index': site_index,
                'scan_number': scan_number,
                'cell_type': cell_type,
                'status': 'success',
                **output_files
            })
        except Exception as e:
            print(f"  ERROR: {e}")
            results.append({
                'site_index': site_index,
                'scan_number': scan_number,
                'cell_type': cell_type,
                'status': f'error: {e}'
            })

    # Summary
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)

    success_count = sum(1 for r in results if r['status'] == 'success')
    print(f"\nSuccessfully generated: {success_count}/{len(SELECTED_SPECTRA)} spectra")

    print("\nOutput files:")
    for r in results:
        if r['status'] == 'success':
            print(f"  {r['site_index']}_{r['scan_number']}:")
            print(f"    PDF:  {os.path.basename(r['pdf'])}")
            print(f"    TIFF: {os.path.basename(r['tiff'])}")
            if r['emf']:
                print(f"    EMF:  {os.path.basename(r['emf'])}")
            else:
                print(f"    SVG:  {os.path.basename(r['svg'])} (EMF conversion failed)")

    print("\n" + "=" * 70)
    print("Done!")
    print("=" * 70)

    return results


if __name__ == "__main__":
    main()
