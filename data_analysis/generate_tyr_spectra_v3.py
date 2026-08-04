#!/usr/bin/env python3
"""
Generate Tyr O-GlcNAc annotated spectra v3.
Updates: z+H formula, no MS2 deisotoping, no isotope matching,
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
# Configuration
# =============================================================================
DATA_BASE = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version'
OUTPUT_BASE = f'{DATA_BASE}/Figures/Tyr_OGlcNAc_spectra_v3'
TOLERANCE_PPM = 20.0
SN_THRESHOLD = 0.0

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
                name = 'HexNAc' if abs(mass - 203.0794) < 0.1 or abs(mass - 299.123) < 0.1 or abs(mass - 528.286) < 0.1 else res + pos
                mods.append({'position': int(pos), 'residue': res, 'mass': mass, 'name': name})
    return mods


def convert_to_protein_positions(best_pos_str, protein_start):
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
    t0 = time.time()

    # Load all Tyr PSMs from ranked files
    base = f'{DATA_BASE}/data_source/point_to_point_response'
    all_tyr = []
    for cell in ['HEK293T', 'HepG2', 'Jurkat']:
        df = pd.read_csv(f'{base}/OGlcNAc_Level1_{cell}_ranked.csv')
        tyr = df[df['modified_residue'] == 'Y'].copy()
        tyr['cell_type'] = cell
        all_tyr.append(tyr)

    tyr_df = pd.concat(all_tyr, ignore_index=True)
    print(f'Total Tyr PSMs: {len(tyr_df)}')
    for cell in ['HEK293T', 'HepG2', 'Jurkat']:
        n = (tyr_df['cell_type'] == cell).sum()
        act = tyr_df[tyr_df['cell_type'] == cell]['Activation_Type'].value_counts().to_dict()
        print(f'  {cell}: {n} ({act})')

    # Convert best positions to protein-level
    tyr_df['Protein_Positions'] = tyr_df.apply(
        lambda row: convert_to_protein_positions(
            str(row.get('Best Positions', '')), row.get('Protein.Start', 1)
        ), axis=1
    )

    # Parse raw file and scan from Spectrum column
    tyr_df['raw_file'] = tyr_df['Spectrum'].apply(lambda s: s.rsplit('.', 3)[0])
    tyr_df['scan_num'] = tyr_df['Spectrum'].apply(lambda s: int(s.rsplit('.', 3)[1]))

    # Create output dirs
    for cell in ['HEK293T', 'HepG2', 'Jurkat']:
        os.makedirs(os.path.join(OUTPUT_BASE, cell), exist_ok=True)

    # Group by raw_file + cell_type so we open each mzML once
    groups = tyr_df.groupby(['raw_file', 'cell_type'])
    n_files = len(groups)
    total_annotated = 0

    for i, ((raw_file, cell_type), psm_group) in enumerate(groups):
        cal_dir = CELL_MZML_DIRS[cell_type]
        cal_mzml = os.path.join(cal_dir, f'{raw_file}_calibrated.mzML')

        if not os.path.exists(cal_mzml):
            # Try without _mz suffix
            alt = os.path.join(cal_dir, f'{raw_file.replace("_mz", "")}_calibrated.mzML')
            if os.path.exists(alt):
                cal_mzml = alt
            else:
                print(f'  SKIP: {cal_mzml} not found')
                continue

        n_psms = len(psm_group)
        print(f'[{i+1}/{n_files}] {raw_file} ({cell_type}): {n_psms} PSMs...', end=' ', flush=True)

        reader = mzml_utils.MzMLReader(cal_mzml)
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
            activation = psm.get('Activation_Type', 'HCD')
            conf = str(psm.get('Confidence.Level', ''))
            gene = str(psm.get('Gene', ''))
            if pd.isna(gene) or gene == 'nan':
                gene = str(psm.get('Protein.ID', 'unknown'))

            prot_pos = str(psm.get('Protein_Positions', ''))
            if pd.isna(prot_pos) or prot_pos == 'nan':
                prot_pos = ''

            site_id = f'{gene}_{prot_pos.replace(";", "_")}' if prot_pos else gene
            safe_site = prot_pos.replace(';', '_') if prot_pos else ''
            out_path = os.path.join(OUTPUT_BASE, cell_type,
                                     f'{gene}_{safe_site}_s{scan}_{conf}_{activation}.pdf')

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

    print(f'\nComplete: {total_annotated}/{len(tyr_df)} spectra in {time.time()-t0:.0f}s')
    print(f'Output: {OUTPUT_BASE}/')


if __name__ == '__main__':
    main()
