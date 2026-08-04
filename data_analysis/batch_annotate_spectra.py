#!/usr/bin/env python3
"""
Batch annotate all OGlcNAc and OGalNAc Level1/Level1b spectra
using calibrated mzML files and GlycoSpectrumAnnotator.

Groups PSMs by raw file to read each calibrated mzML only once.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import re
import time
import traceback
from pathlib import Path
from pyteomics import mzml

from spectrum_annotator_ddzby import SpectrumAnnotator

# =============================================================================
# Configuration
# =============================================================================
CAL_MZML_BASE = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version'
OUTPUT_BASE = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures'
MASTER_LIST = '/tmp/annotation_master_list.csv'
TOLERANCE_PPM = 20.0
SN_THRESHOLD = 0.0  # No S/N filter (manual validation mode)

CELL_MZML_DIRS = {
    'HEK293T': f'{CAL_MZML_BASE}/OGlycoTM_HEK293T',
    'HepG2': f'{CAL_MZML_BASE}/OGlycoTM_HepG2',
    'Jurkat': f'{CAL_MZML_BASE}/OGlycoTM_Jurkat',
}


def parse_modifications(assigned_mods):
    """Parse OPair Assigned.Modifications string."""
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
            pos = ''
            res = ''
            for char in position_part:
                if char.isdigit():
                    pos += char
                else:
                    res += char
            if pos and res:
                mods.append({'position': int(pos), 'residue': res, 'mass': mass, 'name': f'{res}{pos}'})
    return mods


def load_ethcd_scans():
    """Load EThcD scan numbers from the EThcD ranked files.

    Any scan in these files is EThcD; everything else is HCD.
    """
    ethcd_scans = set()
    base = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/point_to_point_response'
    for cell in ['HEK293T', 'HepG2', 'Jurkat']:
        path = os.path.join(base, f'OGlcNAc_Level1_{cell}_EThcD_ranked.csv')
        if os.path.exists(path):
            df = pd.read_csv(path)
            for spec in df['Spectrum']:
                scan = int(spec.split('.')[-3])
                raw_file = spec.rsplit('.', 3)[0]
                ethcd_scans.add((raw_file, scan))
    return ethcd_scans


ETHCD_SCANS = None  # lazy load


def annotate_file(raw_file, cell_type, psm_group, output_base):
    """Annotate all PSMs from one raw file."""
    cal_dir = CELL_MZML_DIRS[cell_type]
    cal_mzml = os.path.join(cal_dir, f'{raw_file}_calibrated.mzML')

    if not os.path.exists(cal_mzml):
        print(f'  WARNING: Calibrated mzML not found: {cal_mzml}')
        return 0

    scan_nums = set(psm_group['scan_num'].tolist())

    # Read spectrum data from calibrated mzML
    spectra = {}
    reader = mzml.MzML(cal_mzml)
    for spec in reader:
        sid = spec.get('id', '')
        if 'scan=' not in sid:
            continue
        scan = int(sid.split('scan=')[-1])
        if scan in scan_nums:
            # Look up activation from pre-built cache
            act_type = ACTIVATION_CACHE.get((raw_file, scan), 'HCD')

            spectra[scan] = {
                'mz': spec['m/z array'],
                'intensity': spec['intensity array'],
                'ms_level': spec.get('ms level', 2),
                'activation': act_type,
            }
            if len(spectra) == len(scan_nums):
                break
    reader.close()

    # Annotate each PSM
    n_annotated = 0
    for _, psm in psm_group.iterrows():
        scan = psm['scan_num']
        if scan not in spectra:
            continue

        spec_data = spectra[scan]
        if len(spec_data['mz']) == 0:
            continue

        peptide = psm['Peptide']
        charge = int(psm['Charge'])
        obs_mz = float(psm['Observed.M.Z'])
        mods = parse_modifications(psm['Assigned.Modifications'])
        glyco_type = psm['glyco_type']

        # Determine gene name
        gene = str(psm.get('Gene', ''))
        if pd.isna(gene) or gene == 'nan':
            gene = str(psm.get('Protein.ID', 'unknown'))

        # Determine activation type from cache
        activation = spec_data.get('activation', 'HCD')

        # Site info
        conf = psm.get('Confidence.Level', '')
        best_pos = str(psm.get('Best Positions', ''))
        if pd.isna(best_pos) or best_pos == 'nan':
            best_pos = ''

        # Gene with site (e.g., "PRDX6 Y89")
        gene_site = f'{gene} {best_pos}' if best_pos else gene

        # Title: Gene Site ActivationType scan# (shown in PDF header)
        title_gene = f'{gene_site} {activation} s{scan}'

        # Output path: Gene_Site_s#_Conf_Activation.pdf
        safe_site = best_pos.replace(';', '_') if best_pos else ''
        out_dir = os.path.join(output_base, f'{glyco_type}_spectra_v2', cell_type)
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, f'{gene}_{safe_site}_s{scan}_{conf}_{activation}.pdf')

        # Skip if already exists
        if os.path.exists(out_path):
            n_annotated += 1
            continue

        try:
            annotator = SpectrumAnnotator(
                peptide=peptide,
                modifications=mods,
                precursor_charge=charge,
                precursor_mz=obs_mz,
                exp_mz=spec_data['mz'],
                exp_intensity=spec_data['intensity'],
                tolerance_ppm=TOLERANCE_PPM,
                site_index=f'{cell_type}_{gene_site}_s{scan}_{conf}',
                gene=title_gene,
                activation_type=activation,
                glycan_labels={},
                do_deisotope=True,
                sn_threshold=SN_THRESHOLD,
            )
            fig = annotator.plot(output_path=out_path)
            plt.close(fig)
            n_annotated += 1
        except Exception as e:
            print(f'    ERROR annotating {gene} s{scan}: {e}')
            traceback.print_exc()

    return n_annotated


ACTIVATION_CACHE = {}  # (raw_file, scan) -> 'HCD' or 'EThcD'
FILTER_CACHE_PATH = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/filter_string_cache.csv'


def main():
    global ACTIVATION_CACHE

    # Load activation type cache
    if os.path.exists(FILTER_CACHE_PATH):
        cache_df = pd.read_csv(FILTER_CACHE_PATH)
        ACTIVATION_CACHE = {(r['raw_file'], int(r['scan'])): r['activation']
                           for _, r in cache_df.iterrows()}
        n_ethcd = sum(1 for v in ACTIVATION_CACHE.values() if v == 'EThcD')
        print(f'Loaded activation cache: {len(ACTIVATION_CACHE)} scans ({n_ethcd} EThcD)')
    else:
        print(f'WARNING: No activation cache at {FILTER_CACHE_PATH}. All scans default to HCD.')

    print(f'Loading master PSM list...')
    master = pd.read_csv(MASTER_LIST)
    print(f'Total PSMs: {len(master)}')

    # Group by raw file and cell type
    groups = master.groupby(['raw_file', 'cell_type'])
    n_files = len(groups)
    print(f'Files to process: {n_files}')

    total_annotated = 0
    total_start = time.time()

    for i, ((raw_file, cell_type), psm_group) in enumerate(groups):
        file_start = time.time()
        n_psms = len(psm_group)
        print(f'[{i+1}/{n_files}] {raw_file} ({cell_type}): {n_psms} PSMs...', end=' ', flush=True)

        try:
            n = annotate_file(raw_file, cell_type, psm_group, OUTPUT_BASE)
            elapsed = time.time() - file_start
            print(f'{n} annotated in {elapsed:.0f}s')
            total_annotated += n
        except Exception as e:
            print(f'FAILED: {e}')
            traceback.print_exc()

    total_elapsed = time.time() - total_start
    print(f'\nDone! {total_annotated} spectra annotated in {total_elapsed/3600:.1f} hours')


if __name__ == '__main__':
    main()
