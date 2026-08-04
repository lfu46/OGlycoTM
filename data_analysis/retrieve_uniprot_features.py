#!/usr/bin/env python3
"""
Retrieve UniProt positional features for O-GalNAc proteins.
Extracts: Signal peptide, Transmembrane, Topological domain, Domain, Repeat.
"""
import json
import time
import pandas as pd
import requests
import os

DATA_BASE = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source'
OUTPUT_JSON = f'{DATA_BASE}/reference/uniprot_features_OGalNAc.json'
OUTPUT_TSV = f'{DATA_BASE}/reference/uniprot_positional_features_OGalNAc.tsv'

FEATURE_TYPES = {
    'Signal', 'Transit peptide', 'Transmembrane',
    'Topological domain', 'Domain', 'Repeat', 'Region',
    'Intramembrane', 'Propeptide', 'Chain',
}


def fetch_uniprot_features(accession, session):
    """Fetch features from UniProt REST API."""
    url = f'https://rest.uniprot.org/uniprotkb/{accession}.json'
    resp = session.get(url, timeout=30)
    if resp.status_code == 200:
        return resp.json()
    return None


def parse_features(data, accession):
    """Parse positional features from UniProt JSON."""
    rows = []
    if not data or 'features' not in data:
        return rows

    for feat in data['features']:
        ftype = feat.get('type', '')
        if ftype not in FEATURE_TYPES:
            continue

        location = feat.get('location', {})
        start = location.get('start', {}).get('value')
        end = location.get('end', {}).get('value')
        desc = feat.get('description', '')

        if start is not None and end is not None:
            rows.append({
                'Protein.ID': accession,
                'feature_type': ftype,
                'start': int(start),
                'end': int(end),
                'description': desc,
            })

    # Also extract protein length
    seq = data.get('sequence', {})
    length = seq.get('length', 0)
    if length:
        rows.append({
            'Protein.ID': accession,
            'feature_type': 'Protein_length',
            'start': 1,
            'end': int(length),
            'description': '',
        })

    return rows


def main():
    # Collect all unique Protein.IDs
    all_pids = set()
    for cell in ['HEK293T', 'HepG2', 'Jurkat']:
        de = pd.read_csv(f'{DATA_BASE}/differential_analysis/OGalNAc_site_DE_{cell}.csv')
        all_pids.update(de['Protein.ID'].unique())
        # Also include protein-level (may have proteins without site DE)
        filt = pd.read_csv(f'{DATA_BASE}/filtered/OGalNAc_{cell}.csv')
        all_pids.update(filt['Protein.ID'].unique())

    all_pids = sorted(all_pids)
    print(f'Total unique proteins: {len(all_pids)}')

    # Check for cached results
    raw_cache = {}
    if os.path.exists(OUTPUT_JSON):
        with open(OUTPUT_JSON) as f:
            raw_cache = json.load(f)
        print(f'Loaded cache: {len(raw_cache)} proteins')

    # Fetch missing
    session = requests.Session()
    session.headers['User-Agent'] = 'OGlycoTM_analysis/1.0 (longping.fu@gatech.edu)'

    to_fetch = [p for p in all_pids if p not in raw_cache]
    print(f'To fetch: {len(to_fetch)} proteins')

    for i, pid in enumerate(to_fetch):
        try:
            data = fetch_uniprot_features(pid, session)
            if data:
                raw_cache[pid] = data
            else:
                raw_cache[pid] = None
                print(f'  {pid}: not found')
        except Exception as e:
            raw_cache[pid] = None
            print(f'  {pid}: error {e}')

        if (i + 1) % 50 == 0:
            print(f'  Fetched {i+1}/{len(to_fetch)}')
            # Save intermediate
            with open(OUTPUT_JSON, 'w') as f:
                json.dump(raw_cache, f)
        time.sleep(0.1)  # rate limit

    # Save final cache
    with open(OUTPUT_JSON, 'w') as f:
        json.dump(raw_cache, f)
    print(f'Saved raw JSON: {len(raw_cache)} proteins')

    # Parse features
    all_features = []
    for pid in all_pids:
        data = raw_cache.get(pid)
        if data:
            rows = parse_features(data, pid)
            all_features.extend(rows)

    feat_df = pd.DataFrame(all_features)
    feat_df.to_csv(OUTPUT_TSV, sep='\t', index=False)
    print(f'Saved parsed features: {len(feat_df)} rows')

    # Summary
    print(f'\nFeature type counts:')
    print(feat_df['feature_type'].value_counts().to_string())

    # Coverage
    has_topo = set(feat_df[feat_df['feature_type'] == 'Topological domain']['Protein.ID'])
    has_tm = set(feat_df[feat_df['feature_type'] == 'Transmembrane']['Protein.ID'])
    has_signal = set(feat_df[feat_df['feature_type'] == 'Signal peptide']['Protein.ID'])
    has_domain = set(feat_df[feat_df['feature_type'] == 'Domain']['Protein.ID'])

    print(f'\nAnnotation coverage ({len(all_pids)} proteins):')
    print(f'  Topological domain: {len(has_topo)}')
    print(f'  Transmembrane: {len(has_tm)}')
    print(f'  Signal peptide: {len(has_signal)}')
    print(f'  Domain: {len(has_domain)}')
    print(f'  Any topology (TOPO_DOM or TM or SP): {len(has_topo | has_tm | has_signal)}')


if __name__ == '__main__':
    main()
