#!/usr/bin/env python3
"""
Annotate all EThcD-supported secretory O-GalNAc sites in Jurkat.
For each site: MS1 isolation + HCD + EThcD for every EThcD PSM.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import re

import mzml_utils
from spectrum_annotator_ddzby import SpectrumAnnotator

DATA_BASE = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version'
OUTPUT_BASE = f'{DATA_BASE}/Figures/OGalNAc_secretory_EThcD_spectra'
TOLERANCE_PPM = 20.0
SN_THRESHOLD = 0.0
FILTER_CACHE_PATH = f'{DATA_BASE}/data_source/filter_string_cache.csv'
MZML_BASE = f'{DATA_BASE}/OGlycoTM_Jurkat'


def parse_modifications(assigned_mods):
    if pd.isna(assigned_mods) or not assigned_mods:
        return []
    mods = []
    for mod in assigned_mods.split(','):
        mod = mod.strip()
        if not mod or '(' not in mod:
            continue
        mass_str = mod[mod.find('(')+1:mod.find(')')]
        try:
            mass = float(mass_str)
        except ValueError:
            continue
        position_part = mod[:mod.find('(')]
        if position_part == 'N-term':
            mods.append({'position': 0, 'residue': 'N-term', 'mass': mass, 'name': 'N-term'})
        elif position_part == 'C-term':
            mods.append({'position': -1, 'residue': 'C-term', 'mass': mass, 'name': 'C-term'})
        else:
            pos, res = '', ''
            for char in position_part:
                if char.isdigit():
                    pos += char
                else:
                    res += char
            if pos and res:
                name = 'HexNAc' if abs(mass - 555.297) < 0.1 or abs(mass - 326.134) < 0.1 or abs(mass - 203.079) < 0.1 or abs(mass - 528.286) < 0.1 or abs(mass - 299.123) < 0.1 else res + pos
                mods.append({'position': int(pos), 'residue': res, 'mass': mass, 'name': name})
    return mods


def get_site_id(row):
    sp = str(row['Site.Probabilities'])
    m = re.match(r'\[(\d+)', sp)
    if not m:
        return ''
    pep_pos = int(m.group(1))
    prot_pos = int(row['Protein.Start']) + pep_pos - 1
    residue = row['Peptide'][pep_pos-1] if pep_pos <= len(row['Peptide']) else '?'
    return f'{row["Gene"]}_{residue}{prot_pos}'


def find_preceding_ms1(reader, scan_num, max_gap=20):
    for s in range(scan_num - 1, max(1, scan_num - max_gap), -1):
        try:
            spec = reader.get_spectrum(s)
            if spec.filter_string and 'Full ms ' in spec.filter_string and 'ms2' not in spec.filter_string.lower():
                return s
        except Exception:
            continue
    return None


def find_paired_hcd(reader, scan_num, obs_mz, max_gap=5, tolerance=0.02):
    for s in range(scan_num - max_gap, scan_num + max_gap + 1):
        if s == scan_num:
            continue
        try:
            spec = reader.get_spectrum(s)
            if spec.precursor_mz and abs(spec.precursor_mz - obs_mz) < tolerance:
                fs = spec.filter_string or ''
                if 'hcd' in fs.lower() and 'etd' not in fs.lower():
                    return s
        except Exception:
            continue
    return None


def annotate_ms1(reader, ms1_scan, precursor_mz, charge, site_id, out_path):
    spec = reader.get_spectrum(ms1_scan)
    if spec.n_peaks == 0:
        return False
    window = 3.0
    mask = (spec.mz >= precursor_mz - window) & (spec.mz <= precursor_mz + window)
    mz_window = spec.mz[mask]
    int_window = spec.intensity[mask]
    if len(mz_window) == 0:
        return False

    tol_da = precursor_mz * 20.0 / 1e6
    close_mask = np.abs(spec.mz - precursor_mz) < tol_da
    if np.any(close_mask):
        closest_idx = np.argmin(np.abs(spec.mz[close_mask] - precursor_mz))
        matched_mz = spec.mz[close_mask][closest_idx]
        mono_ppm = (precursor_mz - matched_mz) / matched_mz * 1e6
    else:
        # Try next MS1
        matched_mz = precursor_mz
        mono_ppm = float('nan')

    fig, ax = plt.subplots(figsize=(8, 3))
    ax.bar(mz_window, int_window, width=0.015, color='grey', alpha=0.8)
    for i in range(min(5, charge + 2)):
        iso_mz = precursor_mz + i * 1.00335 / charge
        ax.axvline(iso_mz, color='red', alpha=0.5, linestyle='--', linewidth=0.8)
    title = f'{site_id} | MS1 scan {ms1_scan} | m/z {precursor_mz:.4f} ({charge}+)'
    title += f' | nearest peak {matched_mz:.4f} ({mono_ppm:+.1f} ppm)'
    ax.set_title(title, fontsize=9)
    ax.set_xlabel('m/z', fontsize=9)
    ax.set_ylabel('Intensity', fontsize=9)
    ax.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return True


def annotate_ms2(reader, scan_num, psm, activation, site_id, out_path):
    spec = reader.get_spectrum(scan_num)
    if spec.n_peaks == 0:
        return False
    mods = parse_modifications(psm['Assigned.Modifications'])
    annotator = SpectrumAnnotator(
        peptide=psm['Peptide'], modifications=mods,
        precursor_charge=int(psm['Charge']), precursor_mz=float(psm['Observed.M.Z']),
        exp_mz=spec.mz, exp_intensity=spec.intensity,
        tolerance_ppm=TOLERANCE_PPM, site_index=site_id, gene=psm['Gene'],
        activation_type=activation, do_deisotope=False,
        sn_threshold=SN_THRESHOLD, scan_num=scan_num,
        confidence_level=str(psm.get('Confidence.Level', '')),
    )
    fig = annotator.plot(output_path=out_path)
    plt.close(fig)
    return True


def main():
    # Load activation cache
    cache_df = pd.read_csv(FILTER_CACHE_PATH)
    act_cache = {(r['raw_file'], int(r['scan'])): r['activation']
                 for _, r in cache_df.iterrows()}
    print(f'Loaded activation cache: {len(act_cache)} scans')

    # Load secretory classification
    features = pd.read_csv(f'{DATA_BASE}/data_source/reference/uniprot_positional_features_OGalNAc.tsv', sep='\t')
    sp_prots = set(features[features['feature_type'] == 'Signal']['Protein.ID'].unique())
    tm_prots = set(features[features['feature_type'] == 'Transmembrane']['Protein.ID'].unique())
    secretory = sp_prots | tm_prots

    # Load PSMs
    psm_df = pd.read_csv(f'{DATA_BASE}/data_source/site/OGalNAc_site_Jurkat.csv')
    psm_df['raw_file'] = psm_df['Spectrum'].apply(lambda s: s.rsplit('.', 3)[0])
    psm_df['scan_num'] = psm_df['Spectrum'].apply(lambda s: int(s.rsplit('.', 3)[1]))
    psm_df['activation'] = psm_df.apply(lambda r: act_cache.get((r['raw_file'], r['scan_num']), 'HCD'), axis=1)
    psm_df['site_id'] = psm_df.apply(get_site_id, axis=1)
    print(f'Loaded {len(psm_df)} PSMs')

    # Filter: EThcD + secretory + Level1/1b
    ethcd_sec = psm_df[
        (psm_df['activation'].str.contains('EThcD')) &
        (psm_df['Protein.ID'].isin(secretory)) &
        (psm_df['Confidence.Level'].isin(['Level1', 'Level1b']))
    ].copy()
    print(f'EThcD secretory PSMs: {len(ethcd_sec)}')
    print(f'Unique sites: {ethcd_sec["site_id"].nunique()}')

    os.makedirs(OUTPUT_BASE, exist_ok=True)
    readers = {}

    for site_id in sorted(ethcd_sec['site_id'].unique()):
        site_psms = ethcd_sec[ethcd_sec['site_id'] == site_id].sort_values('O.Pair.Score', ascending=False)
        print(f'\n=== {site_id} ({len(site_psms)} EThcD PSMs) ===')

        for _, psm in site_psms.iterrows():
            raw_file = psm['raw_file']
            scan = psm['scan_num']
            charge = int(psm['Charge'])
            obs_mz = float(psm['Observed.M.Z'])
            conf = psm['Confidence.Level']

            # Create folder per PSM
            folder = os.path.join(OUTPUT_BASE, f'{site_id}_s{scan}')
            os.makedirs(folder, exist_ok=True)

            # Get reader
            if raw_file not in readers:
                cal_mzml = os.path.join(MZML_BASE, f'{raw_file}_calibrated.mzML')
                if not os.path.exists(cal_mzml):
                    print(f'    SKIP: no mzML for {raw_file}')
                    continue
                readers[raw_file] = mzml_utils.MzMLReader(cal_mzml)
            reader = readers[raw_file]

            # 1. MS1 isolation
            ms1_scan = find_preceding_ms1(reader, scan)
            if ms1_scan is None:
                # Try next MS1
                for s in range(scan + 1, scan + 20):
                    try:
                        spec = reader.get_spectrum(s)
                        if spec.filter_string and 'Full ms ' in spec.filter_string and 'ms2' not in spec.filter_string.lower():
                            ms1_scan = s
                            break
                    except Exception:
                        continue

            if ms1_scan:
                ms1_path = os.path.join(folder, f'{site_id}_s{ms1_scan}_MS1_isolation.pdf')
                if not os.path.exists(ms1_path):
                    # Check if precursor is in this MS1; if not, try next
                    ms1_spec = reader.get_spectrum(ms1_scan)
                    tol_da = obs_mz * 20.0 / 1e6
                    if not np.any(np.abs(ms1_spec.mz - obs_mz) < tol_da):
                        # Try next MS1
                        for s in range(scan + 1, scan + 20):
                            try:
                                sp = reader.get_spectrum(s)
                                if sp.filter_string and 'Full ms ' in sp.filter_string and 'ms2' not in sp.filter_string.lower():
                                    if np.any(np.abs(sp.mz - obs_mz) < tol_da):
                                        ms1_scan = s
                                        ms1_path = os.path.join(folder, f'{site_id}_s{ms1_scan}_MS1_isolation.pdf')
                                        break
                            except Exception:
                                continue
                    annotate_ms1(reader, ms1_scan, obs_mz, charge, site_id, ms1_path)

            # 2. EThcD annotation
            ethcd_path = os.path.join(folder, f'{site_id}_s{scan}_{conf}_EThcD.pdf')
            if not os.path.exists(ethcd_path):
                try:
                    annotate_ms2(reader, scan, psm, 'EThcD', site_id, ethcd_path)
                except Exception as e:
                    print(f'    s{scan} EThcD: ERROR {e}')

            # 3. Paired HCD
            hcd_scan = find_paired_hcd(reader, scan, obs_mz)
            if hcd_scan:
                hcd_path = os.path.join(folder, f'{site_id}_s{hcd_scan}_{conf}_HCD.pdf')
                if not os.path.exists(hcd_path):
                    try:
                        annotate_ms2(reader, hcd_scan, psm, 'HCD', site_id, hcd_path)
                    except Exception as e:
                        print(f'    s{hcd_scan} HCD: ERROR {e}')

            n_files = len([f for f in os.listdir(folder) if f.endswith('.pdf')])
            print(f'    {site_id}_s{scan}/ -> {n_files} files')

    print(f'\nDone! Output: {OUTPUT_BASE}')


if __name__ == '__main__':
    main()
