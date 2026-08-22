#!/usr/bin/env python3
"""
Generate annotated spectra for O-GalNAc Level1/1b PSMs.
Uses updated annotator: z+H formula, no MS2 deisotoping, no isotope matching,
scan/confidence/activation in title.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import re
import time
import traceback

import mzml_utils
from spectrum_annotator_ddzby import SpectrumAnnotator

# =============================================================================
DATA_BASE = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version'
OUTPUT_BASE = f'{DATA_BASE}/Figures/OGalNAc_Level1_spectra'
TOLERANCE_PPM = 20.0
SN_THRESHOLD = 0.0
FILTER_CACHE_PATH = f'{DATA_BASE}/data_source/filter_string_cache.csv'

CELL_MZML_DIRS = {
    'HEK293T': f'{DATA_BASE}/OGlycoTM_HEK293T',
    'HepG2': f'{DATA_BASE}/OGlycoTM_HepG2',
    'Jurkat': f'{DATA_BASE}/OGlycoTM_Jurkat',
}


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
    """Parse site from Site.Probabilities and convert to protein position."""
    if pd.isna(site_prob):
        return '', ''
    m = re.match(r'\[(\d+)', str(site_prob))
    if not m:
        return '', ''
    pep_pos = int(m.group(1))
    prot_pos = int(protein_start) + pep_pos - 1 if pd.notna(protein_start) else pep_pos
    return pep_pos, prot_pos


def main():
    t0 = time.time()

    # Load activation cache
    activation_cache = {}
    if os.path.exists(FILTER_CACHE_PATH):
        cache_df = pd.read_csv(FILTER_CACHE_PATH)
        activation_cache = {(r['raw_file'], int(r['scan'])): r['activation']
                           for _, r in cache_df.iterrows()}
        print(f'Loaded activation cache: {len(activation_cache)} scans')

    # Load OGalNAc Level1/1b PSMs
    all_psms = []
    for cell in ['HEK293T', 'HepG2', 'Jurkat']:
        df = pd.read_csv(f'{DATA_BASE}/data_source/filtered/OGalNAc_{cell}.csv')
        level1 = df[df['Confidence.Level'].isin(['Level1', 'Level1b'])].copy()
        level1['cell_type'] = cell
        all_psms.append(level1)

    psm_df = pd.concat(all_psms, ignore_index=True)
    print(f'Total OGalNAc Level1/1b PSMs: {len(psm_df)}')
    for cell in ['HEK293T', 'HepG2', 'Jurkat']:
        print(f'  {cell}: {(psm_df["cell_type"] == cell).sum()}')

    # Parse raw file and scan
    psm_df['raw_file'] = psm_df['Spectrum'].apply(lambda s: s.rsplit('.', 3)[0])
    psm_df['scan_num'] = psm_df['Spectrum'].apply(lambda s: int(s.rsplit('.', 3)[1]))

    # Parse site positions
    psm_df[['pep_site', 'prot_site']] = psm_df.apply(
        lambda r: pd.Series(parse_site_position(r['Site.Probabilities'], r['Protein.Start'])),
        axis=1
    )

    # Create output dirs
    for cell in ['HEK293T', 'HepG2', 'Jurkat']:
        os.makedirs(os.path.join(OUTPUT_BASE, cell), exist_ok=True)

    # Group by raw_file + cell_type
    groups = psm_df.groupby(['raw_file', 'cell_type'])
    n_files = len(groups)
    total_annotated = 0

    for i, ((raw_file, cell_type), psm_group) in enumerate(groups):
        cal_dir = CELL_MZML_DIRS[cell_type]

        # Try different mzML naming patterns
        cal_mzml = None
        for suffix in ['_calibrated.mzML', '_mz_calibrated.mzML']:
            path = os.path.join(cal_dir, f'{raw_file}{suffix}')
            if os.path.exists(path):
                cal_mzml = path
                break
        if not cal_mzml:
            # Try removing _mz
            alt = os.path.join(cal_dir, f'{raw_file.replace("_mz", "")}_calibrated.mzML')
            if os.path.exists(alt):
                cal_mzml = alt

        if not cal_mzml:
            print(f'  SKIP: no mzML for {raw_file}')
            continue

        n_psms = len(psm_group)
        print(f'[{i+1}/{n_files}] {raw_file} ({cell_type}): {n_psms} PSMs...', end=' ', flush=True)

        reader = mzml_utils.open_spectra(cal_mzml)
        n_ok = 0

        for _, psm in psm_group.iterrows():
            scan = psm['scan_num']
            try:
                spec = reader.get_spectrum(scan)
            except Exception:
                continue
            if spec.n_peaks == 0:
                continue

            peptide = psm['Peptide']
            charge = int(psm['Charge'])
            obs_mz = float(psm['Observed.M.Z'])
            mods = parse_modifications(psm['Assigned.Modifications'])
            conf = str(psm.get('Confidence.Level', ''))

            # Activation type from cache or default HCD
            activation = activation_cache.get((raw_file, scan), 'HCD')

            gene = str(psm.get('Gene', ''))
            if pd.isna(gene) or gene == 'nan':
                gene = str(psm.get('Protein.ID', 'unknown'))

            prot_site = psm.get('prot_site', '')
            pep_site = psm.get('pep_site', '')
            if prot_site and pep_site:
                residue = peptide[int(pep_site)-1] if int(pep_site) <= len(peptide) else '?'
                site_str = f'{residue}{prot_site}'
            else:
                site_str = ''

            site_id = f'{gene}_{site_str}' if site_str else gene
            out_path = os.path.join(OUTPUT_BASE, cell_type,
                                    f'{gene}_{site_str}_s{scan}_{conf}_{activation}.pdf')

            try:
                annotator = SpectrumAnnotator(
                    peptide=peptide,
                    modifications=mods,
                    precursor_charge=charge,
                    precursor_mz=obs_mz,
                    exp_mz=spec.mz,
                    exp_intensity=spec.intensity,
                    tolerance_ppm=TOLERANCE_PPM,
                    site_index=site_id,
                    gene=gene,
                    activation_type=activation,
                    do_deisotope=False,
                    sn_threshold=SN_THRESHOLD,
                    scan_num=scan,
                    confidence_level=conf,
                )
                fig = annotator.plot(output_path=out_path)
                plt.close(fig)
                n_ok += 1
            except Exception as e:
                print(f'\n    ERROR {gene} s{scan}: {e}')
                traceback.print_exc()

        total_annotated += n_ok
        elapsed = time.time() - t0
        print(f'{n_ok}/{n_psms} done ({elapsed:.0f}s total)')

    print(f'\nComplete: {total_annotated}/{len(psm_df)} spectra in {time.time()-t0:.0f}s')
    print(f'Output: {OUTPUT_BASE}/')


if __name__ == '__main__':
    main()
