#!/usr/bin/env python3
"""
Batch annotate non-secretory O-GalNAc EThcD spectra (HEK293T Level1/1b).
Generates MS1 isolation + HCD + EThcD PDFs for each PSM.
Reports MS1 ppm for quality assessment.
"""
import json
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import traceback

from pyteomics import mzml as pyteomics_mzml
from mzml_utils import PROTON, NEUTRON_MASS
from spectrum_annotator_ddzby import SpectrumAnnotator

DATA_BASE = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version'
SOURCE = f'{DATA_BASE}/data_source'
MZML_DIR = f'{DATA_BASE}/OGlycoTM_HEK293T'
OUTPUT_BASE = f'{DATA_BASE}/Figures/OGalNAc_nonsecretory_check'
TOLERANCE_PPM = 20.0
SN_THRESHOLD = 0.0
ISOLATION_WIDTH = 1.4


def parse_modifications(mod_str):
    mods = []
    if pd.isna(mod_str) or not mod_str:
        return mods
    for part in mod_str.split(','):
        part = part.strip()
        if not part or '(' not in part:
            continue
        mass = float(part[part.find('(') + 1:part.find(')')])
        prefix = part[:part.find('(')]
        if prefix == 'N-term':
            mods.append({'position': 0, 'residue': 'N-term', 'mass': mass})
        else:
            pos_str, res = '', ''
            for ch in prefix:
                if ch.isdigit():
                    pos_str += ch
                else:
                    res += ch
            if pos_str and res:
                mods.append({'position': int(pos_str), 'residue': res, 'mass': mass})
    return mods


def find_preceding_ms1(reader, target_scan, precursor_mz=None, max_lookback=500):
    fallback = None
    for sn in range(target_scan - 1, max(1, target_scan - max_lookback), -1):
        try:
            spec = reader.get_spectrum(sn)
            if spec is None or spec.ms_level != 1:
                continue
        except:
            continue
        if fallback is None:
            fallback = spec
        if precursor_mz is None:
            return spec
        tol = precursor_mz * 10e-6
        if ((spec.mz >= precursor_mz - tol) & (spec.mz <= precursor_mz + tol)).any():
            return spec
    return fallback


def plot_ms1(ms1_spec, precursor_mz, charge, neutral_mass, output_path, gene, site_idx, peptide):
    fig, ax = plt.subplots(figsize=(14, 6))
    half_iso = ISOLATION_WIDTH / 2
    margin = 4.0
    mask = (ms1_spec.mz >= precursor_mz - margin) & (ms1_spec.mz <= precursor_mz + margin)
    rmz, rint = ms1_spec.mz[mask], ms1_spec.intensity[mask]
    if len(rmz) == 0:
        plt.close(fig)
        return None
    mx = max(rint)
    rel = rint / mx * 100
    ax.axvspan(precursor_mz - half_iso, precursor_mz + half_iso, alpha=0.08, color='royalblue')
    for e in [precursor_mz - half_iso, precursor_mz + half_iso]:
        ax.axvline(e, color='royalblue', lw=1, ls='--', alpha=0.5)
    ax.vlines(rmz, 0, rel, colors='#b0b0b0', lw=0.8, zorder=1)
    n_iso = min(8, int(neutral_mass / 800) + 3)
    theo = [precursor_mz + n * NEUTRON_MASS / charge for n in range(n_iso)]
    matched = []
    for n, imz in enumerate(theo):
        idx = np.argmin(np.abs(ms1_spec.mz - imz))
        if abs(ms1_spec.mz[idx] - imz) < 0.02:
            emz = ms1_spec.mz[idx]
            ppm = (emz - imz) / imz * 1e6
            erel = ms1_spec.intensity[idx] / mx * 100
            matched.append((n, emz, erel, ppm))
            ax.vlines(emz, 0, erel, colors='#0d75bc', lw=2.5, zorder=3)
            ax.annotate(f'M+{n}' if n else 'M (mono)', (emz, erel),
                        textcoords="offset points", xytext=(0, 8), ha='center',
                        fontsize=8, color='#0d75bc', fontweight='bold')
    ax.set_xlabel('m/z', fontsize=12)
    ax.set_ylabel('Relative Abundance (%)', fontsize=12)
    ax.set_xlim(precursor_mz - margin, precursor_mz + margin)
    ax.set_ylim(0, max(rel) * 1.15)
    ax.set_title(f'{gene} - {site_idx} | MS1 scan {ms1_spec.scan_num} (RT {ms1_spec.rt:.2f} min)\n'
                 f'{peptide} | {precursor_mz:.4f} m/z, z={charge} | M={neutral_mass:.4f} Da', fontsize=10)
    lines = []
    if matched:
        lines.append(f'Isotopes: M+0..M+{max(n for n, _, _, _ in matched)} ({matched[0][3]:.1f} ppm)')
    else:
        lines.append('WARNING: No isotope match')
    lines.append(f'Isolation: {ISOLATION_WIDTH:.1f} m/z')
    ax.text(0.02, 0.97, '\n'.join(lines), transform=ax.transAxes, fontsize=9, va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    return {'ms1_scan': ms1_spec.scan_num, 'mono_ppm': matched[0][3] if matched else np.nan,
            'n_isotopes': len(matched)}


def find_hcd_scan(reader, ethcd_scan):
    ref = reader.get_spectrum(ethcd_scan)
    for offset in [2, 1, 3, 4]:
        candidate = ethcd_scan - offset
        try:
            s = reader.get_spectrum(candidate)
            fs = (s.filter_string or '').lower()
            if 'hcd' in fs and 'etd' not in fs and abs(s.precursor_mz - ref.precursor_mz) < 1.0:
                return candidate
        except:
            pass
    return None


def main():
    os.makedirs(OUTPUT_BASE, exist_ok=True)

    # Load non-secretory classification
    with open(f'{SOURCE}/reference/uniprot_features_OGalNAc.json') as f:
        cache = json.load(f)

    def is_secretory(acc):
        data = cache.get(acc)
        if not data or 'features' not in data:
            return None
        has_sp = any(f['type'] == 'Signal' for f in data['features'])
        has_tm = any(f['type'] == 'Transmembrane' for f in data['features'])
        return has_sp or has_tm

    # Load site data and filter
    site = pd.read_csv(f'{SOURCE}/site/OGalNAc_site_HEK293T.csv')
    site['is_sec'] = site['Protein.ID'].apply(lambda p: is_secretory(p))
    site['scan'] = site['Spectrum'].apply(lambda s: int(s.rsplit('.', 3)[1]))
    site['raw_file'] = site['Spectrum'].apply(lambda s: s.rsplit('.', 3)[0])

    fs_cache = pd.read_csv(f'{SOURCE}/filter_string_cache.csv')
    merged = site.merge(fs_cache[['raw_file', 'scan', 'activation']], on=['raw_file', 'scan'], how='left')

    target = merged[(merged['is_sec'] == False) &
                    (merged['Confidence.Level'].isin(['Level1', 'Level1b'])) &
                    (merged['activation'] == 'EThcD')].copy()

    print(f'PSMs to annotate: {len(target)}')

    summary_rows = []

    for raw_file, group in target.groupby('raw_file'):
        mzml_path = os.path.join(MZML_DIR, f'{raw_file}_calibrated.mzML')
        if not os.path.exists(mzml_path):
            print(f'  SKIP: {mzml_path} not found')
            continue

        # Build set of all scans we need: EThcD + nearby HCD + preceding MS1
        ethcd_scans = set(group['scan'].tolist())
        # Collect a window around each EThcD for HCD and MS1
        need_scans = set()
        for s in ethcd_scans:
            for offset in range(-20, 5):  # MS1 up to 20 before, HCD up to 4 before
                need_scans.add(s + offset)

        # Stream through mzML once, collect needed spectra
        print(f'Streaming {os.path.basename(mzml_path)} ({len(group)} PSMs)...', flush=True)
        spectra = {}
        reader = pyteomics_mzml.MzML(mzml_path)
        for spec in reader:
            sid = spec.get('id', '')
            if 'scan=' not in sid:
                continue
            scan = int(sid.split('scan=')[-1])
            if scan in need_scans:
                mz_array = spec.get('m/z array', np.array([]))
                int_array = spec.get('intensity array', np.array([]))
                ms_level = spec.get('ms level', 0)
                fs = spec.get('filter string', spec.get('scanList', {}).get('scan', [{}])[0].get('filter string', ''))

                precursor_mz = None
                if 'precursorList' in spec:
                    precs = spec['precursorList'].get('precursor', [])
                    if precs:
                        ions = precs[0].get('selectedIonList', {}).get('selectedIon', [])
                        if ions:
                            precursor_mz = ions[0].get('selected ion m/z')
                elif 'selected precursors' in spec:
                    precs = spec['selected precursors']
                    if precs:
                        precursor_mz = precs[0].get('selected ion m/z')

                rt = spec.get('scanList', {}).get('scan', [{}])[0].get('scan start time', 0)

                spectra[scan] = {
                    'mz': mz_array, 'intensity': int_array,
                    'ms_level': ms_level, 'filter_string': str(fs),
                    'precursor_mz': precursor_mz, 'rt': rt,
                }

            # Stop early if we have all needed scans
            if scan > max(need_scans) + 100:
                break

        print(f'  Collected {len(spectra)} spectra', flush=True)

        for _, psm in group.iterrows():
            ethcd_scan = psm['scan']
            gene = str(psm['Gene'])
            peptide = psm['Peptide']
            charge = int(psm['Charge'])
            obs_mz = float(psm['Observed.M.Z'])
            mods = parse_modifications(psm['Assigned.Modifications'])
            conf = str(psm['Confidence.Level'])
            site_idx = psm['site_index']
            neutral_mass = obs_mz * charge - charge * PROTON

            folder = f'{gene}_{site_idx}_s{ethcd_scan}'
            out_dir = os.path.join(OUTPUT_BASE, folder)
            os.makedirs(out_dir, exist_ok=True)

            row = {'site_index': site_idx, 'gene': gene, 'peptide': peptide,
                   'charge': charge, 'obs_mz': obs_mz, 'conf': conf,
                   'ethcd_scan': ethcd_scan, 'raw_file': raw_file,
                   'glycan': psm['Total.Glycan.Composition']}

            # Find paired HCD (scan-1 to scan-4, HCD not EThcD)
            hcd_scan = None
            for offset in [2, 1, 3, 4]:
                candidate = ethcd_scan - offset
                if candidate in spectra:
                    fs = spectra[candidate]['filter_string'].lower()
                    if 'hcd' in fs and 'etd' not in fs:
                        hcd_scan = candidate
                        break

            # Find preceding MS1
            ms1_scan = None
            ref_scan = hcd_scan or ethcd_scan
            for s in range(ref_scan - 1, ref_scan - 20, -1):
                if s in spectra and spectra[s]['ms_level'] == 1:
                    ms1_scan = s
                    break

            # MS1 isolation plot
            if ms1_scan and len(spectra[ms1_scan]['mz']) > 0:
                class MS1Spec:
                    pass
                ms1 = MS1Spec()
                ms1.mz = spectra[ms1_scan]['mz']
                ms1.intensity = spectra[ms1_scan]['intensity']
                ms1.scan_num = ms1_scan
                ms1.rt = spectra[ms1_scan]['rt']

                ms1_out = os.path.join(out_dir, f'{site_idx}_s{ms1_scan}_MS1_isolation.pdf')
                ms1_info = plot_ms1(ms1, obs_mz, charge, neutral_mass, ms1_out,
                                    gene, site_idx, peptide)
                if ms1_info:
                    row['ms1_scan'] = ms1_info['ms1_scan']
                    row['ms1_mono_ppm'] = round(ms1_info['mono_ppm'], 2)
                    row['ms1_n_isotopes'] = ms1_info['n_isotopes']

            # EThcD annotation
            if ethcd_scan in spectra and len(spectra[ethcd_scan]['mz']) > 0:
                ethcd_out = os.path.join(out_dir, f'{site_idx}_s{ethcd_scan}_{conf}_EThcD.pdf')
                try:
                    ann = SpectrumAnnotator(
                        peptide=peptide, modifications=mods, precursor_charge=charge,
                        precursor_mz=obs_mz, exp_mz=spectra[ethcd_scan]['mz'],
                        exp_intensity=spectra[ethcd_scan]['intensity'],
                        tolerance_ppm=TOLERANCE_PPM, site_index=site_idx, gene=gene,
                        activation_type='EThcD', do_deisotope=False, sn_threshold=SN_THRESHOLD,
                        scan_num=ethcd_scan, confidence_level=conf)
                    fig = ann.plot(output_path=ethcd_out)
                    plt.close(fig)
                    fmr = ann.false_match_rate
                    row['ethcd_fmr'] = round(fmr.fmr_peaks, 4)
                    row['ethcd_coverage'] = ann.annotation_stats.get('sequence_coverage_bonds', '')
                except Exception as e:
                    row['ethcd_fmr'] = ''
                    row['ethcd_coverage'] = f'ERROR: {e}'

            # HCD annotation
            if hcd_scan and hcd_scan in spectra and len(spectra[hcd_scan]['mz']) > 0:
                hcd_out = os.path.join(out_dir, f'{site_idx}_s{hcd_scan}_{conf}_HCD.pdf')
                try:
                    ann = SpectrumAnnotator(
                        peptide=peptide, modifications=mods, precursor_charge=charge,
                        precursor_mz=obs_mz, exp_mz=spectra[hcd_scan]['mz'],
                        exp_intensity=spectra[hcd_scan]['intensity'],
                        tolerance_ppm=TOLERANCE_PPM, site_index=site_idx, gene=gene,
                        activation_type='HCD', do_deisotope=False, sn_threshold=SN_THRESHOLD,
                        scan_num=hcd_scan, confidence_level=conf)
                    fig = ann.plot(output_path=hcd_out)
                    plt.close(fig)
                    fmr = ann.false_match_rate
                    row['hcd_scan'] = hcd_scan
                    row['hcd_fmr'] = round(fmr.fmr_peaks, 4)
                    row['hcd_coverage'] = ann.annotation_stats.get('sequence_coverage_bonds', '')
                except Exception as e:
                    row['hcd_fmr'] = ''
                    row['hcd_coverage'] = f'ERROR: {e}'

            summary_rows.append(row)

        print(f'  {raw_file}: {len(group)} PSMs done', flush=True)

    # Save summary
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(os.path.join(OUTPUT_BASE, 'annotation_summary.csv'), index=False)

    print(f'\nDone: {len(summary_rows)} PSMs annotated')
    print(f'Output: {OUTPUT_BASE}/')

    # MS1 ppm summary
    if 'ms1_mono_ppm' in summary.columns:
        ppm = summary['ms1_mono_ppm'].dropna()
        print(f'\nMS1 mono ppm: median={ppm.median():.1f}, mean={ppm.mean():.1f}, '
              f'range=[{ppm.min():.1f}, {ppm.max():.1f}]')
        print(f'  |ppm| < 5: {(ppm.abs() < 5).sum()}/{len(ppm)}')
        print(f'  |ppm| < 10: {(ppm.abs() < 10).sum()}/{len(ppm)}')


if __name__ == '__main__':
    main()
