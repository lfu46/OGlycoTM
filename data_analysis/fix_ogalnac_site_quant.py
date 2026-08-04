#!/usr/bin/env python3
"""
Fix O-GalNAc site-level quantification for EThcD scans.
EThcD scans have weak TMT reporters — use paired HCD scan reporters instead.

For each EThcD PSM in the site file:
  1. Find the paired HCD scan (same precursor, adjacent scan number)
  2. If HCD PSM exists in OPair output → use its TMT intensities
  3. Otherwise → extract TMT reporters from the HCD mzML spectrum directly
"""
import pandas as pd
import numpy as np
import os
import mzml_utils

DATA_BASE = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version'
CACHE_PATH = f'{DATA_BASE}/data_source/filter_string_cache.csv'
CELL_MZML_DIRS = {
    'HEK293T': f'{DATA_BASE}/OGlycoTM_HEK293T',
    'HepG2': f'{DATA_BASE}/OGlycoTM_HepG2',
    'Jurkat': f'{DATA_BASE}/OGlycoTM_Jurkat',
}

# TMT6plex reporter m/z values
TMT_REPORTERS = {
    'Intensity.Tuni_1': 126.127726,
    'Intensity.Tuni_2': 127.131081,
    'Intensity.Tuni_3': 128.134436,
    'Intensity.Ctrl_4': 129.137790,
    'Intensity.Ctrl_5': 130.141145,
    'Intensity.Ctrl_6': 131.138180,
}
TMT_TOL_DA = 0.01  # 10 mDa tolerance for reporter extraction

INTENSITY_COLS = list(TMT_REPORTERS.keys())


def extract_tmt_reporters(reader, scan_num):
    """Extract TMT reporter ion intensities from an mzML spectrum."""
    try:
        spec = reader.get_spectrum(scan_num)
    except Exception:
        return None

    mz, intensity = spec.mz, spec.intensity
    if len(mz) == 0:
        return None

    reporters = {}
    for col, theo_mz in TMT_REPORTERS.items():
        mask = np.abs(mz - theo_mz) < TMT_TOL_DA
        if mask.any():
            reporters[col] = float(intensity[mask].max())
        else:
            reporters[col] = 0.0
    return reporters


def find_paired_hcd_scan(reader, ethcd_scan, precursor_mz):
    """Find the HCD scan paired with an EThcD scan."""
    for offset in [-1, -2, -3, -4]:
        hcd_scan = ethcd_scan + offset
        try:
            spec = reader.get_spectrum(hcd_scan)
            fs = str(spec.filter_string) if hasattr(spec, 'filter_string') and spec.filter_string else ''
            if 'etd' in fs.lower():
                continue  # Skip other EThcD scans
            if 'Full ms ' in fs and 'ms2' not in fs:
                continue  # Skip MS1 scans
            if hasattr(spec, 'precursor_mz') and spec.precursor_mz:
                if abs(spec.precursor_mz - precursor_mz) < 0.01:
                    return hcd_scan
        except Exception:
            continue
    return None


def main():
    # Load activation cache
    cache = pd.read_csv(CACHE_PATH)
    act_dict = {(r['raw_file'], int(r['scan'])): r['activation'] for _, r in cache.iterrows()}

    for cell in ['HEK293T', 'HepG2', 'Jurkat']:
        print(f"\n{'='*60}")
        print(f"  {cell}")
        print(f"{'='*60}")

        site_path = f'{DATA_BASE}/data_source/site/OGalNAc_site_{cell}.csv'
        filtered_path = f'{DATA_BASE}/data_source/filtered/OGalNAc_{cell}.csv'

        site_df = pd.read_csv(site_path)
        filtered_df = pd.read_csv(filtered_path)

        # Parse raw_file and scan
        site_df['raw_file'] = site_df['Spectrum'].apply(lambda s: s.rsplit('.', 3)[0])
        site_df['scan_num'] = site_df['Spectrum'].apply(lambda s: int(s.rsplit('.', 3)[1]))
        site_df['activation'] = site_df.apply(
            lambda r: act_dict.get((r['raw_file'], r['scan_num']), 'HCD'), axis=1
        )

        # Build OPair lookup from filtered data
        filtered_df['raw_file'] = filtered_df['Spectrum'].apply(lambda s: s.rsplit('.', 3)[0])
        filtered_df['scan_num'] = filtered_df['Spectrum'].apply(lambda s: int(s.rsplit('.', 3)[1]))
        opair_lookup = {}
        for _, r in filtered_df.iterrows():
            opair_lookup[(r['raw_file'], r['scan_num'])] = r

        ethcd_mask = site_df['activation'] == 'EThcD'
        n_ethcd = ethcd_mask.sum()
        n_total = len(site_df)
        print(f"  Site PSMs: {n_total} total, {n_ethcd} EThcD")

        if n_ethcd == 0:
            print("  No EThcD PSMs to fix.")
            continue

        # Process EThcD PSMs grouped by raw file
        readers = {}
        n_from_opair = 0
        n_from_mzml = 0
        n_failed = 0

        for idx in site_df[ethcd_mask].index:
            row = site_df.loc[idx]
            raw = row['raw_file']
            scan = row['scan_num']
            prec_mz = row['Observed.M.Z']

            # Try to find paired HCD in OPair output first
            paired_hcd = None
            for offset in [-1, -2, -3, -4]:
                hcd_scan = scan + offset
                key = (raw, hcd_scan)
                if key in opair_lookup:
                    hcd_row = opair_lookup[key]
                    if abs(hcd_row['Observed.M.Z'] - prec_mz) < 0.01:
                        # Use OPair HCD intensities
                        for col in INTENSITY_COLS:
                            site_df.at[idx, col] = hcd_row[col]
                        paired_hcd = hcd_scan
                        n_from_opair += 1
                        break

            if paired_hcd is None:
                # Read from mzML directly
                if raw not in readers:
                    for suffix in ['_calibrated.mzML', '_mz_calibrated.mzML']:
                        path = os.path.join(CELL_MZML_DIRS[cell], f'{raw}{suffix}')
                        if os.path.exists(path):
                            readers[raw] = mzml_utils.MzMLReader(path)
                            break

                if raw in readers:
                    hcd_scan = find_paired_hcd_scan(readers[raw], scan, prec_mz)
                    if hcd_scan:
                        reporters = extract_tmt_reporters(readers[raw], hcd_scan)
                        if reporters:
                            for col in INTENSITY_COLS:
                                site_df.at[idx, col] = reporters[col]
                            n_from_mzml += 1
                        else:
                            n_failed += 1
                    else:
                        n_failed += 1
                else:
                    n_failed += 1

        print(f"  Fixed: {n_from_opair} from OPair HCD, {n_from_mzml} from mzML, {n_failed} failed")

        # Save updated site file (backup original first)
        backup_path = site_path.replace('.csv', '_original.csv')
        if not os.path.exists(backup_path):
            import shutil
            shutil.copy(site_path, backup_path)
            print(f"  Backed up original: {backup_path}")

        # Drop temp columns before saving
        site_df.drop(columns=['raw_file', 'scan_num', 'activation'], inplace=True)
        site_df.to_csv(site_path, index=False)
        print(f"  Saved: {site_path}")

    print("\nDone!")


if __name__ == '__main__':
    main()
