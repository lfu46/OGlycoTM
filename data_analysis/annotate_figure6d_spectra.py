#!/usr/bin/env python3
"""
Generate MS1 + HCD + EThcD annotated spectra for Figure 6D selected sites.
Each site gets its own subfolder with all available spectrum types.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import re
import traceback

import mzml_utils
from spectrum_annotator_ddzby import SpectrumAnnotator

DATA_BASE = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version'
OUTPUT_BASE = f'{DATA_BASE}/Figures/Figure6D_secretory_spectra_review'
TOLERANCE_PPM = 20.0
SN_THRESHOLD = 0.0
FILTER_CACHE_PATH = f'{DATA_BASE}/data_source/filter_string_cache.csv'
MZML_DIR = f'{DATA_BASE}/OGlycoTM_Jurkat'

# Target sites: (gene, residue+number, preferred_scan or None for best)
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


def parse_site_position(site_prob, protein_start):
    if pd.isna(site_prob):
        return '', ''
    m = re.match(r'\[(\d+)', str(site_prob))
    if not m:
        return '', ''
    pep_pos = int(m.group(1))
    prot_pos = int(protein_start) + pep_pos - 1 if pd.notna(protein_start) else pep_pos
    return pep_pos, prot_pos


def find_preceding_ms1(reader, scan_num, max_gap=20):
    """Find the MS1 scan preceding the given MS2 scan."""
    for s in range(scan_num - 1, max(1, scan_num - max_gap), -1):
        try:
            spec = reader.get_spectrum(s)
            if spec.filter_string and 'Full ms ' in spec.filter_string and 'ms2' not in spec.filter_string.lower():
                return s
        except Exception:
            continue
    return None


def find_paired_ethcd(reader, activation_cache, raw_file, hcd_scan, precursor_mz, max_gap=5, tolerance=0.02):
    """Find an EThcD scan paired with the given HCD scan (same precursor)."""
    for s in range(hcd_scan - max_gap, hcd_scan + max_gap + 1):
        if s == hcd_scan:
            continue
        act = activation_cache.get((raw_file, s), '')
        if 'EThcD' in act:
            try:
                spec = reader.get_spectrum(s)
                if spec.precursor_mz and abs(spec.precursor_mz - precursor_mz) < tolerance:
                    return s
            except Exception:
                continue
    return None


def annotate_ms1(reader, ms1_scan, precursor_mz, charge, site_id, out_path):
    """Generate MS1 isolation window plot."""
    spec = reader.get_spectrum(ms1_scan)
    if spec.n_peaks == 0:
        return False

    # Plot MS1 around precursor
    window = 3.0  # Da
    mask = (spec.mz >= precursor_mz - window) & (spec.mz <= precursor_mz + window)
    mz_window = spec.mz[mask]
    int_window = spec.intensity[mask]

    if len(mz_window) == 0:
        return False

    # Find closest peak to precursor m/z for ppm calculation
    tol_da = precursor_mz * 20.0 / 1e6  # 20 ppm window
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

    # Highlight theoretical isotope peaks based on observed precursor m/z
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
    """Generate annotated MS2 spectrum."""
    spec = reader.get_spectrum(scan_num)
    if spec.n_peaks == 0:
        return False

    peptide = psm['Peptide']
    charge = int(psm['Charge'])
    obs_mz = float(psm['Observed.M.Z'])
    mods = parse_modifications(psm['Assigned.Modifications'])
    conf = str(psm.get('Confidence.Level', ''))

    annotator = SpectrumAnnotator(
        peptide=peptide,
        modifications=mods,
        precursor_charge=charge,
        precursor_mz=obs_mz,
        exp_mz=spec.mz,
        exp_intensity=spec.intensity,
        tolerance_ppm=TOLERANCE_PPM,
        site_index=site_id,
        gene=psm['Gene'],
        activation_type=activation,
        do_deisotope=False,
        sn_threshold=SN_THRESHOLD,
        scan_num=scan_num,
        confidence_level=conf,
    )
    fig = annotator.plot(output_path=out_path)
    plt.close(fig)
    return True


def main():
    # Load activation cache
    activation_cache = {}
    if os.path.exists(FILTER_CACHE_PATH):
        cache_df = pd.read_csv(FILTER_CACHE_PATH)
        activation_cache = {(r['raw_file'], int(r['scan'])): r['activation']
                           for _, r in cache_df.iterrows()}
        print(f'Loaded activation cache: {len(activation_cache)} scans')

    # Load OGalNAc Jurkat PSMs
    psm_df = pd.read_csv(f'{DATA_BASE}/data_source/site/OGalNAc_site_Jurkat.csv')
    psm_df['raw_file'] = psm_df['Spectrum'].apply(lambda s: s.rsplit('.', 3)[0])
    psm_df['scan_num'] = psm_df['Spectrum'].apply(lambda s: int(s.rsplit('.', 3)[1]))
    psm_df[['pep_site', 'prot_site']] = psm_df.apply(
        lambda r: pd.Series(parse_site_position(r['Site.Probabilities'], r['Protein.Start'])),
        axis=1
    )
    print(f'Loaded {len(psm_df)} Jurkat OGalNAc PSMs')

    # Open readers cache
    readers = {}

    for gene, site_str in TARGETS:
        residue = site_str[0]
        number = int(site_str[1:])
        site_id = f'{gene}_{site_str}'
        print(f'\n=== {site_id} ===')

        # Find best PSM for this site
        gene_psms = psm_df[psm_df['Gene'] == gene].copy()
        site_psms = gene_psms[gene_psms['Site.Probabilities'].str.contains(f'{residue}{number}', na=False)]
        if len(site_psms) == 0:
            site_psms = gene_psms

        # Prefer Level1 with highest OPair score
        level1 = site_psms[site_psms['Confidence.Level'].isin(['Level1', 'Level1b'])]
        if len(level1) > 0:
            # Prefer EThcD scans if available
            ethcd_psms = []
            for _, p in level1.iterrows():
                act = activation_cache.get((p['raw_file'], p['scan_num']), 'HCD')
                if 'EThcD' in act:
                    ethcd_psms.append(p)

            if ethcd_psms:
                best = max(ethcd_psms, key=lambda p: p['O.Pair.Score'])
            else:
                best = level1.sort_values('O.Pair.Score', ascending=False).iloc[0]
        else:
            best = site_psms.sort_values('O.Pair.Score', ascending=False).iloc[0]

        raw_file = best['raw_file']
        scan = best['scan_num']
        charge = int(best['Charge'])
        obs_mz = float(best['Observed.M.Z'])
        activation = activation_cache.get((raw_file, scan), 'HCD')

        print(f'  Best PSM: {raw_file} scan {scan} ({activation})')
        print(f'  Peptide: {best["Peptide"]}  charge: {charge}  m/z: {obs_mz:.4f}')
        print(f'  Level: {best["Confidence.Level"]}  OPair: {best["O.Pair.Score"]:.3f}')

        # Create output folder
        folder = os.path.join(OUTPUT_BASE, f'{site_id}_s{scan}')
        os.makedirs(folder, exist_ok=True)

        # Get reader
        if raw_file not in readers:
            cal_mzml = os.path.join(MZML_DIR, f'{raw_file}_calibrated.mzML')
            if not os.path.exists(cal_mzml):
                print(f'  SKIP: no mzML for {raw_file}')
                continue
            readers[raw_file] = mzml_utils.MzMLReader(cal_mzml)
        reader = readers[raw_file]

        # 1. MS1 isolation check
        ms1_scan = find_preceding_ms1(reader, scan)
        if ms1_scan:
            ms1_path = os.path.join(folder, f'{site_id}_s{ms1_scan}_MS1_isolation.pdf')
            ok = annotate_ms1(reader, ms1_scan, obs_mz, charge, site_id, ms1_path)
            print(f'  MS1 scan {ms1_scan}: {"OK" if ok else "FAILED"}')
        else:
            print(f'  MS1: not found')

        # 2. Annotate the best scan (HCD or EThcD)
        ms2_path = os.path.join(folder, f'{site_id}_s{scan}_{best["Confidence.Level"]}_{activation}.pdf')
        ok = annotate_ms2(reader, scan, best, activation, site_id, ms2_path)
        print(f'  {activation} scan {scan}: {"OK" if ok else "FAILED"}')

        # 3. If best is HCD, also look for paired EThcD
        if activation == 'HCD':
            ethcd_scan = find_paired_ethcd(reader, activation_cache, raw_file, scan, obs_mz)
            if ethcd_scan:
                ethcd_path = os.path.join(folder, f'{site_id}_s{ethcd_scan}_{best["Confidence.Level"]}_EThcD.pdf')
                ok = annotate_ms2(reader, ethcd_scan, best, 'EThcD', site_id, ethcd_path)
                print(f'  EThcD scan {ethcd_scan}: {"OK" if ok else "FAILED"}')
            else:
                print(f'  EThcD: no paired scan found')

        # 4. If best is EThcD, also look for paired HCD
        if 'EThcD' in activation:
            # Look for HCD in nearby scans
            for s in range(scan - 5, scan + 6):
                if s == scan:
                    continue
                act = activation_cache.get((raw_file, s), '')
                if act == 'HCD':
                    try:
                        spec = reader.get_spectrum(s)
                        if spec.precursor_mz and abs(spec.precursor_mz - obs_mz) < 0.02:
                            hcd_path = os.path.join(folder, f'{site_id}_s{s}_{best["Confidence.Level"]}_HCD.pdf')
                            ok = annotate_ms2(reader, s, best, 'HCD', site_id, hcd_path)
                            print(f'  HCD scan {s}: {"OK" if ok else "FAILED"}')
                            break
                    except Exception:
                        continue

        # 5. Also check if there are other good PSMs with different activation
        other_psms = level1 if len(level1) > 0 else site_psms
        for _, alt_psm in other_psms.iterrows():
            alt_scan = alt_psm['scan_num']
            alt_raw = alt_psm['raw_file']
            if alt_scan == scan and alt_raw == raw_file:
                continue
            alt_act = activation_cache.get((alt_raw, alt_scan), 'HCD')
            # Only annotate if it's a different activation type
            if alt_act != activation:
                if alt_raw not in readers:
                    cal_mzml = os.path.join(MZML_DIR, f'{alt_raw}_calibrated.mzML')
                    if os.path.exists(cal_mzml):
                        readers[alt_raw] = mzml_utils.MzMLReader(cal_mzml)
                if alt_raw in readers:
                    alt_path = os.path.join(folder, f'{site_id}_s{alt_scan}_{alt_psm["Confidence.Level"]}_{alt_act}.pdf')
                    if not os.path.exists(alt_path):
                        ok = annotate_ms2(readers[alt_raw], alt_scan, alt_psm, alt_act, site_id, alt_path)
                        print(f'  Alt {alt_act} scan {alt_scan}: {"OK" if ok else "FAILED"}')
                        break

    print(f'\n=== Done! Output: {OUTPUT_BASE} ===')


if __name__ == '__main__':
    main()
