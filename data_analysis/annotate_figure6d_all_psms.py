#!/usr/bin/env python3
"""
Annotate ALL Level1/1b PSMs for Figure 6D sites.
Each PSM gets its own folder: {site}_{scan}/ with MS1 + HCD + EThcD.
PTPRC capped at top 5 HCD + all EThcD.
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
OUTPUT_BASE = f'{DATA_BASE}/Figures/Figure6D_secretory_spectra_review'
TOLERANCE_PPM = 20.0
SN_THRESHOLD = 0.0
FILTER_CACHE_PATH = f'{DATA_BASE}/data_source/filter_string_cache.csv'
MZML_BASE = f'{DATA_BASE}/OGlycoTM_Jurkat'

TARGETS = [
    ('TFRC', 'T104'),
    ('PTPRC', 'T139'),
    ('PTPRC', 'S146'),
    ('IGSF8', 'T169'),
    ('CPD', 'T44'),
    ('MAN1A2', 'T139'),
    ('B3GALT6', 'S46'),
    ('OS9', 'S529'),
    ('EDEM1', 'S48'),
]


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


def find_preceding_ms1(reader, scan_num, max_gap=20):
    for s in range(scan_num - 1, max(1, scan_num - max_gap), -1):
        try:
            spec = reader.get_spectrum(s)
            if spec.filter_string and 'Full ms ' in spec.filter_string and 'ms2' not in spec.filter_string.lower():
                return s
        except Exception:
            continue
    return None


def find_paired_scan(reader, scan_num, obs_mz, target_act='HCD', max_gap=5, tolerance=0.02):
    for s in range(scan_num - max_gap, scan_num + max_gap + 1):
        if s == scan_num:
            continue
        try:
            spec = reader.get_spectrum(s)
            if spec.precursor_mz and abs(spec.precursor_mz - obs_mz) < tolerance:
                fs = spec.filter_string or ''
                if target_act == 'HCD' and 'hcd' in fs.lower() and 'etd' not in fs.lower():
                    return s
                elif target_act == 'EThcD' and 'etd' in fs.lower():
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
    activation_cache = {}
    if os.path.exists(FILTER_CACHE_PATH):
        cache_df = pd.read_csv(FILTER_CACHE_PATH)
        activation_cache = {(r['raw_file'], int(r['scan'])): r['activation']
                           for _, r in cache_df.iterrows()}
        print(f'Loaded activation cache: {len(activation_cache)} scans')

    psm_df = pd.read_csv(f'{DATA_BASE}/data_source/site/OGalNAc_site_Jurkat.csv')
    psm_df['raw_file'] = psm_df['Spectrum'].apply(lambda s: s.rsplit('.', 3)[0])
    psm_df['scan_num'] = psm_df['Spectrum'].apply(lambda s: int(s.rsplit('.', 3)[1]))
    print(f'Loaded {len(psm_df)} PSMs')

    readers = {}

    for gene, site_str in TARGETS:
        residue = site_str[0]
        number = site_str[1:]
        site_id = f'{gene}_{site_str}'
        print(f'\n=== {site_id} ===')

        gene_psms = psm_df[psm_df['Gene'] == gene]
        site_psms = gene_psms[gene_psms['Site.Probabilities'].str.contains(f'{residue}{number}', na=False)]
        if len(site_psms) == 0:
            site_psms = gene_psms
        level1 = site_psms[site_psms['Confidence.Level'].isin(['Level1', 'Level1b'])].copy()
        level1['activation'] = level1.apply(
            lambda r: activation_cache.get((r['raw_file'], r['scan_num']), 'HCD'), axis=1)

        # Cap PTPRC at top 5 HCD + all EThcD
        if gene == 'PTPRC':
            ethcd = level1[level1['activation'].str.contains('EThcD')]
            hcd = level1[~level1['activation'].str.contains('EThcD')].sort_values('O.Pair.Score', ascending=False).head(5)
            level1 = pd.concat([ethcd, hcd]).drop_duplicates(subset=['scan_num', 'raw_file'])

        print(f'  {len(level1)} PSMs to annotate')

        for _, psm in level1.iterrows():
            raw_file = psm['raw_file']
            scan = psm['scan_num']
            activation = psm['activation']
            conf = psm['Confidence.Level']
            charge = int(psm['Charge'])
            obs_mz = float(psm['Observed.M.Z'])

            # Each PSM gets its own folder
            folder = os.path.join(OUTPUT_BASE, f'{site_id}_s{scan}')
            os.makedirs(folder, exist_ok=True)

            if raw_file not in readers:
                cal_mzml = os.path.join(MZML_BASE, f'{raw_file}_calibrated.mzML')
                if not os.path.exists(cal_mzml):
                    print(f'    SKIP: no mzML for {raw_file}')
                    continue
                readers[raw_file] = mzml_utils.open_spectra(cal_mzml)
            reader = readers[raw_file]

            # 1. MS1 isolation
            ms1_scan = find_preceding_ms1(reader, scan)
            if ms1_scan:
                ms1_path = os.path.join(folder, f'{site_id}_s{ms1_scan}_MS1_isolation.pdf')
                if not os.path.exists(ms1_path):
                    annotate_ms1(reader, ms1_scan, obs_mz, charge, site_id, ms1_path)

            # 2. MS2 (the PSM scan itself)
            ms2_path = os.path.join(folder, f'{site_id}_s{scan}_{conf}_{activation}.pdf')
            if not os.path.exists(ms2_path):
                try:
                    annotate_ms2(reader, scan, psm, activation, site_id, ms2_path)
                except Exception as e:
                    print(f'    s{scan} {activation}: ERROR {e}')

            # 3. Paired scan
            if activation == 'HCD':
                paired = find_paired_scan(reader, scan, obs_mz, 'EThcD')
                if paired:
                    p_path = os.path.join(folder, f'{site_id}_s{paired}_{conf}_EThcD.pdf')
                    if not os.path.exists(p_path):
                        try:
                            annotate_ms2(reader, paired, psm, 'EThcD', site_id, p_path)
                        except Exception:
                            pass
            elif 'EThcD' in activation:
                paired = find_paired_scan(reader, scan, obs_mz, 'HCD')
                if paired:
                    p_path = os.path.join(folder, f'{site_id}_s{paired}_{conf}_HCD.pdf')
                    if not os.path.exists(p_path):
                        try:
                            annotate_ms2(reader, paired, psm, 'HCD', site_id, p_path)
                        except Exception:
                            pass

            n_files = len([f for f in os.listdir(folder) if f.endswith('.pdf')])
            print(f'    {site_id}_s{scan}/ -> {n_files} files')

    # Clean up old loose PDFs in root
    for f in os.listdir(OUTPUT_BASE):
        fp = os.path.join(OUTPUT_BASE, f)
        if f.endswith('.pdf') and os.path.isfile(fp):
            os.remove(fp)

    print(f'\nDone! Output: {OUTPUT_BASE}')


if __name__ == '__main__':
    main()
