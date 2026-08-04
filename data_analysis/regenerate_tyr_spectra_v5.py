#!/usr/bin/env python3
"""
Regenerate Tyr O-GlcNAc annotated spectra v5 for 4 selected sites.
Uses GlycoSpectrumAnnotator v0.3.0 (spectrum_annotator_ddzby).
Generates MS1 isolation + HCD + EThcD annotated PDFs + validation summary.
"""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import traceback

import mzml_utils
from mzml_utils import PROTON, NEUTRON_MASS
from spectrum_annotator_ddzby import SpectrumAnnotator

# =============================================================================
# Configuration
# =============================================================================
DATA_BASE = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version'
EXTRACTED_BASE = f'{DATA_BASE}/data_source/point_to_point_response'
OUTPUT_BASE = f'{DATA_BASE}/Figures/Tyr_OGlcNAc_spectra_v5'

TOLERANCE_PPM = 20.0
SN_THRESHOLD = 0.0  # 0 for manual validation
ISOLATION_WIDTH = 1.4  # Th, typical for Orbitrap

CELL_MZML_DIRS = {
    'HEK293T': f'{DATA_BASE}/OGlycoTM_HEK293T',
    'HepG2': f'{DATA_BASE}/OGlycoTM_HepG2',
    'Jurkat': f'{DATA_BASE}/OGlycoTM_Jurkat',
}

# 4 selected sites with full PSM info from extracted_spectra_EThcD summary
SELECTED_SITES = [
    {
        'name': 'PRDX6_Y89',
        'cell': 'Jurkat',
        'ethcd_scan': 21346,
        'peptide': 'DINAYNCEEPTEK',
        'charge': 3,
        'obs_mz': 761.7037,  # calibrated
        'mods_json': [
            {"position": 0, "residue": "N-term", "mass": 229.1629},
            {"position": 5, "residue": "Y", "mass": 299.123},
            {"position": 13, "residue": "K", "mass": 229.1629},
        ],
        'gene': 'PRDX6',
        'protein_id': 'P30041',
        'site_index': 'P30041_Y89',
        'confidence': 'Level1',
        'spectrum_file': 'P30041_Y89_21346.csv',
        'mzml_file': 'Eclipse_LF_OGlycoTM_Jurkat_OG_15_10222025_calibrated.mzML',
    },
    {
        'name': 'DDX17_Y580',
        'cell': 'HEK293T',
        'ethcd_scan': 22821,
        'peptide': 'TTSSANNPNLMYQDECDRR',
        'charge': 3,
        'obs_mz': 939.4349,  # calibrated
        'mods_json': [
            {"position": 0, "residue": "N-term", "mass": 229.1629},
            {"position": 11, "residue": "M", "mass": 15.9949},
            {"position": 12, "residue": "Y", "mass": 299.123},
            {"position": 16, "residue": "C", "mass": 57.0215},
        ],
        'gene': 'DDX17',
        'protein_id': 'Q92841',
        'site_index': 'Q92841_Y580',
        'confidence': 'Level1',
        'spectrum_file': 'Q92841_Y580_22821.csv',
        'mzml_file': 'Eclipse_LF_OGlycoTM_HEK293T_OG_22_10062025_calibrated.mzML',
    },
    {
        'name': 'HPRT1_Y105',
        'cell': 'HEK293T',
        'ethcd_scan': 25463,
        'peptide': 'SYCNDQSTGDIK',
        'charge': 3,
        'obs_mz': 773.0611,  # calibrated
        'mods_json': [
            {"position": 0, "residue": "N-term", "mass": 229.1629},
            {"position": 2, "residue": "Y", "mass": 528.2859},
            {"position": 12, "residue": "K", "mass": 229.1629},
        ],
        'gene': 'HPRT1',
        'protein_id': 'P00492',
        'site_index': 'P00492_Y105',
        'confidence': 'Level1',
        'spectrum_file': 'P00492_Y105_25463.csv',
        'mzml_file': 'Eclipse_LF_OGlycoTM_HEK293T_OG_5_10062025_mz_calibrated.mzML',
    },
    {
        'name': 'SON_Y270',
        'cell': 'HepG2',
        'ethcd_scan': 43812,
        'peptide': 'SSEPVVTMSVEYQMK',
        'charge': 3,
        'obs_mz': 906.8073,  # calibrated
        'mods_json': [
            {"position": 0, "residue": "N-term", "mass": 229.1629},
            {"position": 8, "residue": "M", "mass": 15.9949},
            {"position": 12, "residue": "Y", "mass": 528.2859},
            {"position": 15, "residue": "K", "mass": 229.1629},
        ],
        'gene': 'SON',
        'protein_id': 'P18583',
        'site_index': 'P18583_Y270',
        'confidence': 'Level1',
        'spectrum_file': 'P18583_Y270_43812.csv',
        'mzml_file': 'Eclipse_LF_OGlycoTM_HepG2_OG_2_11212025_calibrated.mzML',
    },
]


def find_hcd_scan(reader, ethcd_scan):
    """Find the paired HCD scan for an EThcD scan (typically 2 scans before)."""
    for offset in [2, 1, 3, 4]:
        candidate = ethcd_scan - offset
        try:
            spec = reader.get_spectrum(candidate)
            fs = spec.filter_string or ''
            if 'hcd' in fs.lower() and 'etd' not in fs.lower():
                return candidate
        except Exception:
            continue
    return None


def find_preceding_ms1(reader, target_scan, precursor_mz=None, max_lookback=500):
    """Find the MS1 scan preceding the target that contains the precursor."""
    fallback = None
    for scan_num in range(target_scan - 1, max(1, target_scan - max_lookback), -1):
        try:
            spec = reader.get_spectrum(scan_num)
            if spec is None or spec.ms_level != 1:
                continue
        except Exception:
            continue
        if fallback is None:
            fallback = spec
        if precursor_mz is None:
            return spec
        tol = precursor_mz * 10e-6
        mask = (spec.mz >= precursor_mz - tol) & (spec.mz <= precursor_mz + tol)
        if mask.any():
            return spec
    return fallback


def plot_ms1(ms1_spec, precursor_mz, charge, neutral_mass, isolation_width,
             output_path, gene, site_index, peptide, glycan_mass=0.0):
    """MS1 isolation window plot with isotope envelope and co-isolated species."""
    fig, ax = plt.subplots(figsize=(14, 6))
    half_iso = isolation_width / 2
    peptide_mass = neutral_mass - glycan_mass

    margin = 4.0
    mask = ((ms1_spec.mz >= precursor_mz - margin) &
            (ms1_spec.mz <= precursor_mz + margin))
    region_mz = ms1_spec.mz[mask]
    region_int = ms1_spec.intensity[mask]

    if len(region_mz) == 0:
        plt.close(fig)
        return None

    max_int = max(region_int)
    rel_int = region_int / max_int * 100

    # Isolation window shading
    ax.axvspan(precursor_mz - half_iso, precursor_mz + half_iso,
               alpha=0.08, color='royalblue', zorder=0)
    for edge in [precursor_mz - half_iso, precursor_mz + half_iso]:
        ax.axvline(edge, color='royalblue', linewidth=1, linestyle='--', alpha=0.5)

    # All peaks in region (grey)
    ax.vlines(region_mz, 0, rel_int, colors='#b0b0b0', linewidth=0.8, zorder=1)

    # Theoretical isotope envelope
    n_isotopes = min(8, int(neutral_mass / 800) + 3)
    theo_isotopes = [precursor_mz + n * NEUTRON_MASS / charge
                     for n in range(n_isotopes)]

    matched_isotopes = []
    for n, iso_mz in enumerate(theo_isotopes):
        idx = np.argmin(np.abs(ms1_spec.mz - iso_mz))
        if abs(ms1_spec.mz[idx] - iso_mz) < 0.02:
            exp_mz = ms1_spec.mz[idx]
            ppm = (exp_mz - iso_mz) / iso_mz * 1e6
            exp_rel = ms1_spec.intensity[idx] / max_int * 100
            matched_isotopes.append((n, exp_mz, exp_rel, ppm))
            ax.vlines(exp_mz, 0, exp_rel, colors='#0d75bc', linewidth=2.5, zorder=3)
            label = f'M+{n}' if n > 0 else 'M (mono)'
            ax.annotate(label, (exp_mz, exp_rel),
                        textcoords="offset points", xytext=(0, 8),
                        ha='center', fontsize=8, color='#0d75bc', fontweight='bold')

    # Co-isolated species
    iso_mask = ((ms1_spec.mz >= precursor_mz - half_iso) &
                (ms1_spec.mz <= precursor_mz + half_iso))
    coisolated = []
    for mz_val, int_val in zip(ms1_spec.mz[iso_mask], ms1_spec.intensity[iso_mask]):
        is_target = any(abs(mz_val - iso) < 0.02 for iso in theo_isotopes)
        rel = int_val / max_int * 100
        if not is_target and rel > 2.0:
            coisolated.append((mz_val, rel))
            ax.vlines(mz_val, 0, rel, colors='#e74c3c', linewidth=2, zorder=2)
            ax.annotate(f'{mz_val:.2f}', (mz_val, rel),
                        textcoords="offset points", xytext=(0, 8),
                        ha='center', fontsize=7, color='#e74c3c')

    ax.set_xlabel('m/z', fontsize=12)
    ax.set_ylabel('Relative Abundance (%)', fontsize=12)
    ax.set_xlim(precursor_mz - margin, precursor_mz + margin)
    ax.set_ylim(0, max(rel_int) * 1.15)

    title = (f'{gene} - {site_index} | MS1 scan {ms1_spec.scan_num} '
             f'(RT {ms1_spec.rt:.2f} min)\n'
             f'{peptide} | Precursor: {precursor_mz:.4f} m/z, z={charge} | '
             f'M = {neutral_mass:.4f} Da')
    ax.set_title(title, fontsize=10)

    # Summary text
    lines = []
    if matched_isotopes:
        iso_range = f'M+0 to M+{max(n for n, _, _, _ in matched_isotopes)}'
        mono_ppm = matched_isotopes[0][3]
        lines.append(f'Isotopes matched: {iso_range} ({mono_ppm:.1f} ppm mono)')
    else:
        lines.append('WARNING: No isotope peaks matched')
    if coisolated:
        co_str = ', '.join(f'{mz:.2f} ({rel:.0f}%)' for mz, rel in coisolated[:5])
        lines.append(f'Co-isolated ({len(coisolated)}): {co_str}')
    else:
        lines.append('No co-isolated species')
    lines.append(f'Isolation window: {isolation_width:.1f} m/z')
    ax.text(0.02, 0.97, '\n'.join(lines), transform=ax.transAxes, fontsize=9,
            va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))

    handles = [
        mpatches.Patch(color='#0d75bc', alpha=0.8, label='Target isotopes'),
        mpatches.Patch(color='#e74c3c', alpha=0.8, label='Co-isolated'),
        mpatches.Patch(color='royalblue', alpha=0.15, label='Isolation window'),
    ]
    ax.legend(handles=handles, loc='upper right', fontsize=9)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)

    return {
        'ms1_scan': ms1_spec.scan_num,
        'n_isotopes_matched': len(matched_isotopes),
        'mono_ppm': matched_isotopes[0][3] if matched_isotopes else np.nan,
        'n_coisolated': len(coisolated),
    }


def main():
    os.makedirs(OUTPUT_BASE, exist_ok=True)

    readers = {}
    validation_rows = []

    for site in SELECTED_SITES:
        name = site['name']
        cell = site['cell']
        ethcd_scan = site['ethcd_scan']
        charge = site['charge']
        obs_mz = site['obs_mz']
        mzml_path = os.path.join(CELL_MZML_DIRS[cell], site['mzml_file'])
        neutral_mass = obs_mz * charge - charge * PROTON

        print(f'\n=== {name} ({cell}) ===')

        out_dir = os.path.join(OUTPUT_BASE, name)
        os.makedirs(out_dir, exist_ok=True)

        if mzml_path not in readers:
            print(f'  Opening {os.path.basename(mzml_path)}...')
            readers[mzml_path] = mzml_utils.MzMLReader(mzml_path)
        reader = readers[mzml_path]

        row = {
            'site': name,
            'cell_type': cell,
            'gene': site['gene'],
            'peptide': site['peptide'],
            'charge': charge,
            'precursor_mz': obs_mz,
            'neutral_mass': neutral_mass,
            'confidence': site['confidence'],
        }

        # === MS1 isolation window ===
        hcd_scan = find_hcd_scan(reader, ethcd_scan)
        ms1_spec = find_preceding_ms1(reader, hcd_scan or ethcd_scan, obs_mz)
        if ms1_spec is not None:
            ms1_out = os.path.join(out_dir,
                                    f'{name}_s{ms1_spec.scan_num}_MS1_isolation.pdf')
            ms1_info = plot_ms1(ms1_spec, obs_mz, charge, neutral_mass,
                                ISOLATION_WIDTH, ms1_out,
                                site['gene'], site['site_index'],
                                site['peptide'])
            if ms1_info:
                print(f'  MS1 scan {ms1_info["ms1_scan"]}: '
                      f'{ms1_info["n_isotopes_matched"]} isotopes, '
                      f'{ms1_info["mono_ppm"]:.1f} ppm mono, '
                      f'{ms1_info["n_coisolated"]} co-isolated')
                row.update({
                    'ms1_scan': ms1_info['ms1_scan'],
                    'ms1_isotopes_matched': ms1_info['n_isotopes_matched'],
                    'ms1_mono_ppm': round(ms1_info['mono_ppm'], 2),
                    'ms1_n_coisolated': ms1_info['n_coisolated'],
                })
        else:
            print(f'  WARNING: No MS1 scan found')

        # === EThcD annotation ===
        csv_path = os.path.join(EXTRACTED_BASE,
                                f'extracted_spectra_EThcD_{cell}',
                                site['spectrum_file'])
        spec_df = pd.read_csv(csv_path)
        mz = spec_df['mz'].values
        intensity = spec_df['intensity'].values
        print(f'  EThcD scan {ethcd_scan}: {len(mz)} peaks')

        ethcd_out = os.path.join(out_dir,
                                  f'{name}_s{ethcd_scan}_{site["confidence"]}_EThcD.pdf')
        try:
            ethcd_annotator = SpectrumAnnotator(
                peptide=site['peptide'],
                modifications=site['mods_json'],
                precursor_charge=charge,
                precursor_mz=obs_mz,
                exp_mz=mz,
                exp_intensity=intensity,
                tolerance_ppm=TOLERANCE_PPM,
                site_index=site['site_index'],
                gene=site['gene'],
                activation_type='EThcD',
                do_deisotope=False,
                sn_threshold=SN_THRESHOLD,
                scan_num=ethcd_scan,
                confidence_level=site['confidence'],
            )
            fig = ethcd_annotator.plot(output_path=ethcd_out)
            plt.close(fig)
            print(f'  -> {os.path.basename(ethcd_out)}')

            # Validation stats
            ethcd_stats = ethcd_annotator.annotation_stats
            ethcd_fmr = ethcd_annotator.false_match_rate
            cz_ions = [m for m in ethcd_annotator.matched_ions
                       if m.ion_type in ('c', 'z')]
            c_ions = sorted([m for m in cz_ions if m.ion_type == 'c'],
                            key=lambda x: x.ion_number)
            z_ions = sorted([m for m in cz_ions if m.ion_type == 'z'],
                            key=lambda x: x.ion_number)

            row.update({
                'ethcd_scan': ethcd_scan,
                'ethcd_seq_coverage': ethcd_stats.get('sequence_coverage', ''),
                'ethcd_seq_bonds': ethcd_stats.get('sequence_coverage_bonds', ''),
                'ethcd_peaks_annotated': ethcd_stats.get('peaks_annotated', ''),
                'ethcd_intensity_annotated': ethcd_stats.get('intensity_annotated', ''),
                'ethcd_fmr_peaks': round(ethcd_fmr.fmr_peaks, 4) if ethcd_fmr else '',
                'ethcd_fmr_intensity': round(ethcd_fmr.fmr_intensity, 4) if ethcd_fmr else '',
                'ethcd_matched_peaks': ethcd_fmr.matched_peaks if ethcd_fmr else '',
                'ethcd_c_ions': ','.join(
                    f"c{m.ion_number}{'*' if m.has_modification else ''}"
                    for m in c_ions) if c_ions else 'none',
                'ethcd_z_ions': ','.join(
                    f"z{m.ion_number}{'*' if m.has_modification else ''}"
                    for m in z_ions) if z_ions else 'none',
            })
            print(f'    EThcD coverage: {ethcd_stats.get("sequence_coverage_bonds", "?")}, '
                  f'FMR: {ethcd_fmr.fmr_peaks*100:.1f}%')
            print(f'    c-ions: {row["ethcd_c_ions"]}')
            print(f'    z-ions: {row["ethcd_z_ions"]}')
        except Exception as e:
            print(f'  ERROR EThcD: {e}')
            traceback.print_exc()

        # === HCD annotation ===
        if hcd_scan:
            hcd_spec = reader.get_spectrum(hcd_scan)
            print(f'  HCD scan {hcd_scan}: {hcd_spec.n_peaks} peaks')

            hcd_out = os.path.join(out_dir,
                                    f'{name}_s{hcd_scan}_{site["confidence"]}_HCD.pdf')
            try:
                hcd_annotator = SpectrumAnnotator(
                    peptide=site['peptide'],
                    modifications=site['mods_json'],
                    precursor_charge=charge,
                    precursor_mz=obs_mz,
                    exp_mz=hcd_spec.mz,
                    exp_intensity=hcd_spec.intensity,
                    tolerance_ppm=TOLERANCE_PPM,
                    site_index=site['site_index'],
                    gene=site['gene'],
                    activation_type='HCD',
                    do_deisotope=False,
                    sn_threshold=SN_THRESHOLD,
                    scan_num=hcd_scan,
                    confidence_level=site['confidence'],
                )
                fig = hcd_annotator.plot(output_path=hcd_out)
                plt.close(fig)
                print(f'  -> {os.path.basename(hcd_out)}')

                hcd_stats = hcd_annotator.annotation_stats
                hcd_fmr = hcd_annotator.false_match_rate
                row.update({
                    'hcd_scan': hcd_scan,
                    'hcd_seq_coverage': hcd_stats.get('sequence_coverage', ''),
                    'hcd_seq_bonds': hcd_stats.get('sequence_coverage_bonds', ''),
                    'hcd_peaks_annotated': hcd_stats.get('peaks_annotated', ''),
                    'hcd_intensity_annotated': hcd_stats.get('intensity_annotated', ''),
                    'hcd_fmr_peaks': round(hcd_fmr.fmr_peaks, 4) if hcd_fmr else '',
                    'hcd_fmr_intensity': round(hcd_fmr.fmr_intensity, 4) if hcd_fmr else '',
                    'hcd_matched_peaks': hcd_fmr.matched_peaks if hcd_fmr else '',
                })
                print(f'    HCD coverage: {hcd_stats.get("sequence_coverage_bonds", "?")}, '
                      f'FMR: {hcd_fmr.fmr_peaks*100:.1f}%')
            except Exception as e:
                print(f'  ERROR HCD: {e}')
                traceback.print_exc()
        else:
            print(f'  WARNING: No paired HCD scan found for EThcD {ethcd_scan}')

        validation_rows.append(row)

    # Save validation summary
    val_df = pd.DataFrame(validation_rows)
    val_path = os.path.join(OUTPUT_BASE, 'validation_summary.csv')
    val_df.to_csv(val_path, index=False)
    print(f'\nValidation summary: {val_path}')

    # Print summary table
    print(f'\n{"="*80}')
    print('VALIDATION SUMMARY')
    print(f'{"="*80}')
    for r in validation_rows:
        print(f'\n{r["site"]} ({r["cell_type"]})')
        print(f'  Peptide: {r["peptide"]}, z={r["charge"]}, '
              f'M/Z={r["precursor_mz"]:.4f}')
        if 'ms1_mono_ppm' in r:
            print(f'  MS1: {r["ms1_isotopes_matched"]} isotopes, '
                  f'{r["ms1_mono_ppm"]:.1f} ppm, '
                  f'{r["ms1_n_coisolated"]} co-isolated')
        if 'hcd_fmr_peaks' in r:
            print(f'  HCD: coverage {r["hcd_seq_bonds"]}, '
                  f'FMR {float(r["hcd_fmr_peaks"])*100:.1f}%')
        if 'ethcd_fmr_peaks' in r:
            print(f'  EThcD: coverage {r["ethcd_seq_bonds"]}, '
                  f'FMR {float(r["ethcd_fmr_peaks"])*100:.1f}%')
            print(f'    c-ions: {r.get("ethcd_c_ions", "")}')
            print(f'    z-ions: {r.get("ethcd_z_ions", "")}')

    print(f'\nDone. Output: {OUTPUT_BASE}/')


if __name__ == '__main__':
    main()
