#!/usr/bin/env python3
"""
Annotate O-GalNAc sites with glycodomain structural classes.
Uses UniProt positional features to derive per-site:
  - topology_class: extracellular / luminal / cytosolic / secreted / unknown
  - domain_context: inside_domain / outside_domain / repeat_region
  - tm_distance: residues to nearest TM helix
  - has_signal_peptide: boolean
"""
import pandas as pd
import numpy as np
import os

DATA_BASE = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source'
FEATURES_TSV = f'{DATA_BASE}/reference/uniprot_positional_features_OGalNAc.tsv'
OUTPUT = f'{DATA_BASE}/annotation/OGalNAc_site_glycodomain_annotation.csv'


def load_features():
    """Load parsed UniProt features into per-protein lookup."""
    df = pd.read_csv(FEATURES_TSV, sep='\t')
    features = {}
    for pid, group in df.groupby('Protein.ID'):
        features[pid] = group.to_dict('records')
    return features


def get_features_by_type(protein_features, ftype):
    """Get all features of a given type for a protein."""
    return [f for f in protein_features if f['feature_type'] == ftype]


def classify_topology(site_number, protein_features):
    """Determine topology class for a site."""
    topo_doms = get_features_by_type(protein_features, 'Topological domain')
    for td in topo_doms:
        if td['start'] <= site_number <= td['end']:
            desc = str(td['description']).lower()
            if 'extracellular' in desc:
                return 'extracellular'
            elif 'lumen' in desc or 'luminal' in desc:
                return 'luminal'
            elif 'cytoplasmic' in desc or 'cytosolic' in desc:
                return 'cytoplasmic'

    # Check if inside transmembrane
    tms = get_features_by_type(protein_features, 'Transmembrane')
    for tm in tms:
        if tm['start'] <= site_number <= tm['end']:
            return 'transmembrane'

    # No TOPO_DOM match — check if protein is secreted (signal peptide, no TM)
    has_sp = len(get_features_by_type(protein_features, 'Signal')) > 0
    has_tm = len(tms) > 0

    if has_sp and not has_tm:
        return 'secreted'

    # Has TM but site not in any annotated TOPO_DOM — infer from position
    if has_tm and topo_doms:
        # Some TOPO_DOMs may not cover the full sequence; site may be in a gap
        return 'unknown'

    return 'unknown'


def classify_domain_context(site_number, protein_features):
    """Determine whether site is inside a domain, repeat, or linker/stalk."""
    domains = get_features_by_type(protein_features, 'Domain')
    for d in domains:
        if d['start'] <= site_number <= d['end']:
            return 'inside_domain'

    repeats = get_features_by_type(protein_features, 'Repeat')
    for r in repeats:
        if r['start'] <= site_number <= r['end']:
            return 'repeat_region'

    return 'outside_domain'


def compute_tm_distance(site_number, protein_features):
    """Compute minimum distance to nearest TM helix boundary."""
    tms = get_features_by_type(protein_features, 'Transmembrane')
    if not tms:
        return np.nan

    distances = []
    for tm in tms:
        d = min(abs(site_number - tm['start']), abs(site_number - tm['end']))
        # If inside TM, distance is 0
        if tm['start'] <= site_number <= tm['end']:
            d = 0
        distances.append(d)

    return min(distances)


def main():
    os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)

    features = load_features()
    print(f'Loaded features for {len(features)} proteins')

    # Load all site DE data
    all_sites = []
    for cell in ['HEK293T', 'HepG2', 'Jurkat']:
        de = pd.read_csv(f'{DATA_BASE}/differential_analysis/OGalNAc_site_DE_{cell}.csv')
        de['CellType'] = cell
        all_sites.append(de[['site_index', 'Protein.ID', 'Gene', 'site_number',
                             'modified_residue', 'logFC', 'adj.P.Val', 'CellType']])

    sites = pd.concat(all_sites, ignore_index=True)
    print(f'Total sites: {len(sites)}')

    # Annotate each site
    results = []
    for _, row in sites.iterrows():
        pid = row['Protein.ID']
        site_num = int(row['site_number'])
        pf = features.get(pid, [])

        has_sp = len(get_features_by_type(pf, 'Signal')) > 0
        has_tm = len(get_features_by_type(pf, 'Transmembrane')) > 0

        topo = classify_topology(site_num, pf) if pf else 'unknown'
        domain = classify_domain_context(site_num, pf) if pf else 'unknown'
        tm_dist = compute_tm_distance(site_num, pf)

        results.append({
            'site_index': row['site_index'],
            'Protein.ID': pid,
            'Gene': row['Gene'],
            'site_number': site_num,
            'modified_residue': row['modified_residue'],
            'CellType': row['CellType'],
            'logFC': row['logFC'],
            'adj.P.Val': row['adj.P.Val'],
            'topology_class': topo,
            'domain_context': domain,
            'tm_distance': tm_dist,
            'has_signal_peptide': has_sp,
            'has_transmembrane': has_tm,
        })

    result_df = pd.DataFrame(results)
    result_df.to_csv(OUTPUT, index=False)
    print(f'\nSaved: {OUTPUT}')

    # Summary
    print(f'\n=== Topology class ===')
    print(result_df['topology_class'].value_counts().to_string())

    print(f'\n=== Domain context ===')
    print(result_df['domain_context'].value_counts().to_string())

    print(f'\n=== By cell type ===')
    for cell in ['HEK293T', 'HepG2', 'Jurkat']:
        sub = result_df[result_df['CellType'] == cell]
        down = sub[(sub['logFC'] < -0.5) & (sub['adj.P.Val'] < 0.05)]
        print(f'\n{cell}: {len(sub)} sites, {len(down)} downregulated')
        print(f'  Topology: {sub["topology_class"].value_counts().to_dict()}')
        print(f'  Down topology: {down["topology_class"].value_counts().to_dict()}')
        print(f'  Domain: {sub["domain_context"].value_counts().to_dict()}')
        print(f'  Down domain: {down["domain_context"].value_counts().to_dict()}')

    # Spot checks
    print(f'\n=== Spot checks ===')
    for check in [('P02786', 104, 'TFRC T104'),
                   ('P08575', 146, 'PTPRC S146'),
                   ('Q13438', 529, 'OS9 S529')]:
        pid, sn, label = check
        match = result_df[(result_df['Protein.ID'] == pid) & (result_df['site_number'] == sn)]
        if len(match) > 0:
            r = match.iloc[0]
            print(f'  {label}: topology={r["topology_class"]}, domain={r["domain_context"]}, '
                  f'tm_dist={r["tm_distance"]}, SP={r["has_signal_peptide"]}')
        else:
            print(f'  {label}: not in site DE')


if __name__ == '__main__':
    main()
