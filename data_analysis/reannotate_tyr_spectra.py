#!/usr/bin/env python3
"""
Re-annotate Tyrosine O-GlcNAc spectra with corrected protein-level positions.
Fixes the issue where Best Positions used peptide-level positions (e.g., Y12)
instead of protein-level positions (e.g., Y580).

Protein position = Protein.Start + peptide_position - 1
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
# Configuration
# =============================================================================
CAL_MZML_BASE = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version'
OUTPUT_BASE = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/Tyr_OGlcNAc_spectra_v2'
MASTER_LIST = '/tmp/annotation_master_list.csv'
TOLERANCE_PPM = 20.0
SN_THRESHOLD = 0.0

CELL_MZML_DIRS = {
    'HEK293T': f'{CAL_MZML_BASE}/OGlycoTM_HEK293T',
    'HepG2': f'{CAL_MZML_BASE}/OGlycoTM_HepG2',
    'Jurkat': f'{CAL_MZML_BASE}/OGlycoTM_Jurkat',
}

ACTIVATION_CACHE = {}
FILTER_CACHE_PATH = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/filter_string_cache.csv'


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


def convert_to_protein_positions(best_pos_str, protein_start):
    """Convert peptide-level positions to protein-level positions.

    e.g., "Y12" with Protein.Start=569 -> "Y580"
          "T3;S4;Y5" with Protein.Start=86 -> "T88;S89;Y90"
    """
    if not best_pos_str or pd.isna(best_pos_str) or best_pos_str == 'nan':
        return ''

    protein_start = int(protein_start)
    parts = best_pos_str.split(';')
    converted = []
    for part in parts:
        part = part.strip()
        match = re.match(r'([A-Z])(\d+)', part)
        if match:
            residue = match.group(1)
            pep_pos = int(match.group(2))
            prot_pos = protein_start + pep_pos - 1
            converted.append(f'{residue}{prot_pos}')
        else:
            converted.append(part)
    return ';'.join(converted)


def main():
    global ACTIVATION_CACHE

    # Load activation cache
    if os.path.exists(FILTER_CACHE_PATH):
        cache_df = pd.read_csv(FILTER_CACHE_PATH)
        ACTIVATION_CACHE = {(r['raw_file'], int(r['scan'])): r['activation']
                           for _, r in cache_df.iterrows()}
        n_ethcd = sum(1 for v in ACTIVATION_CACHE.values() if v == 'EThcD')
        print(f'Loaded activation cache: {len(ACTIVATION_CACHE)} scans ({n_ethcd} EThcD)')

    # Load master list
    print('Loading master PSM list...')
    master = pd.read_csv(MASTER_LIST)
    print(f'Total PSMs: {len(master)}')

    # Filter for Tyr O-GlcNAc only
    # A PSM has Tyr modification if the peptide has Y at the Best Positions
    tyr_mask = master['Best Positions'].astype(str).str.contains('Y', na=False)
    oglcnac_mask = master['glyco_type'] == 'OGlcNAc'
    localized_mask = master['Confidence.Level'].isin(['Level1', 'Level1b'])

    tyr_psms = master[tyr_mask & oglcnac_mask & localized_mask].copy()
    print(f'Tyr O-GlcNAc PSMs to re-annotate: {len(tyr_psms)}')

    for cell in ['HEK293T', 'HepG2', 'Jurkat']:
        n = (tyr_psms['cell_type'] == cell).sum()
        print(f'  {cell}: {n}')

    # Convert Best Positions to protein positions
    tyr_psms['Protein_Positions'] = tyr_psms.apply(
        lambda row: convert_to_protein_positions(
            str(row.get('Best Positions', '')),
            row.get('Protein.Start', 1)
        ), axis=1
    )

    # Show some examples
    print('\nPosition conversion examples:')
    for _, row in tyr_psms.head(10).iterrows():
        print(f'  {row["Gene"]} {row["Best Positions"]} -> {row["Protein_Positions"]} '
              f'(Protein.Start={row["Protein.Start"]}, peptide={row["Peptide"][:20]}...)')

    # Group by raw file and cell type
    groups = tyr_psms.groupby(['raw_file', 'cell_type'])
    n_files = len(groups)
    print(f'\nFiles to process: {n_files}')

    total_annotated = 0
    total_start = time.time()

    for i, ((raw_file, cell_type), psm_group) in enumerate(groups):
        file_start = time.time()
        cal_dir = CELL_MZML_DIRS[cell_type]
        cal_mzml = os.path.join(cal_dir, f'{raw_file}_calibrated.mzML')

        if not os.path.exists(cal_mzml):
            print(f'  WARNING: Calibrated mzML not found: {cal_mzml}')
            continue

        n_psms = len(psm_group)
        print(f'[{i+1}/{n_files}] {raw_file} ({cell_type}): {n_psms} PSMs...', end=' ', flush=True)

        # Open indexed mzML reader
        reader = mzml_utils.MzMLReader(cal_mzml)

        # Annotate each PSM using indexed access
        n_annotated = 0
        for _, psm in psm_group.iterrows():
            scan = psm['scan_num']
            try:
                spec_obj = reader.get_spectrum(scan)
            except Exception:
                continue

            if spec_obj.n_peaks == 0:
                continue

            activation = ACTIVATION_CACHE.get((raw_file, scan), 'HCD')
            exp_mz = spec_obj.mz
            exp_intensity = spec_obj.intensity

            peptide = psm['Peptide']
            charge = int(psm['Charge'])
            obs_mz = float(psm['Observed.M.Z'])
            mods = parse_modifications(psm['Assigned.Modifications'])

            gene = str(psm.get('Gene', ''))
            if pd.isna(gene) or gene == 'nan':
                gene = str(psm.get('Protein.ID', 'unknown'))

            conf = psm.get('Confidence.Level', '')

            # Use PROTEIN positions (the fix)
            prot_pos = str(psm.get('Protein_Positions', ''))
            if pd.isna(prot_pos) or prot_pos == 'nan':
                prot_pos = ''

            gene_site = f'{gene} {prot_pos}' if prot_pos else gene
            title_gene = f'{gene_site} {activation} s{scan}'

            safe_site = prot_pos.replace(';', '_') if prot_pos else ''
            out_dir = os.path.join(OUTPUT_BASE, cell_type)
            os.makedirs(out_dir, exist_ok=True)
            out_path = os.path.join(out_dir, f'{gene}_{safe_site}_s{scan}_{conf}_{activation}.pdf')

            # Remove old file with wrong position if it exists
            # (old pattern: Gene_Y{peptide_pos}_s{scan}_...)

            try:
                annotator = SpectrumAnnotator(
                    peptide=peptide,
                    modifications=mods,
                    precursor_charge=charge,
                    precursor_mz=obs_mz,
                    exp_mz=exp_mz,
                    exp_intensity=exp_intensity,
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

        reader.close()
        elapsed = time.time() - file_start
        print(f'{n_annotated} annotated in {elapsed:.0f}s')
        total_annotated += n_annotated

    total_elapsed = time.time() - total_start
    print(f'\nDone! {total_annotated} Tyr spectra re-annotated in {total_elapsed/60:.1f} minutes')

    # Summary
    for cell in ['HEK293T', 'HepG2', 'Jurkat']:
        d = os.path.join(OUTPUT_BASE, cell)
        if os.path.exists(d):
            n = len([f for f in os.listdir(d) if f.endswith('.pdf')])
            print(f'  {cell}: {n} files')


if __name__ == '__main__':
    main()
