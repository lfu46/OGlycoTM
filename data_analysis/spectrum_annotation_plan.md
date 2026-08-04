# Spectral Annotation Workflow for O-GlcNAc Site Localization

## Objective
Generate high-quality annotated EThcD spectra demonstrating confident O-GlcNAcylation site localization for manuscript revision. The final selected examples must:
1. Show excellent spectral quality with clear c/z ion coverage
2. Have peptides with multiple S/T sites to demonstrate localization confidence
3. Fit the biological story in Figure 6 (structured vs IDR site responsiveness)
4. Ideally fall within functional protein domains for mechanistic explanation

---

## Phase 1: Data Preparation and Quality Assessment

### Step 1.1: Load and parse O-Pair/FragPipe search results

**Input files needed (user to provide paths):**
- O-Pair PSM output file (typically `*_psm.tsv` from FragPipe)
- Protein database used for search (FASTA)

**Key columns to extract:**
- Peptide sequence
- Protein accession
- Modification site position(s)
- O-Pair confidence level (level 1, 1b, 2, 3)
- O-Pair score
- Hyperscore
- Delta score (difference between best and second-best localization)
- Expectation value (E-value)
- Spectrum reference (scan number, file name)
- Precursor m/z, charge
- Retention time

**Code approach:**
```python
import pandas as pd

# Load PSM file
psm_df = pd.read_csv(psm_file_path, sep='\t')

# Display column names to identify relevant scoring columns
print(psm_df.columns.tolist())

# Basic statistics on scoring distributions
scoring_cols = ['OPair Score', 'Hyperscore', 'Expectation', 'Delta Score']  # adjust names as needed
psm_df[scoring_cols].describe()
```

### Step 1.2: Apply quality filters

**Recommended filtering criteria:**
```python
# Filter 1: O-Pair confidence level (CRITICAL - only level 1 or 1b)
filtered_df = psm_df[psm_df['Localization Level'].isin(['1', '1b', 1, '1.0'])]

# Filter 2: Minimum O-Pair score threshold (explore distribution first)
# Suggest: top 25% of scores as starting point
score_threshold = filtered_df['OPair Score'].quantile(0.75)
filtered_df = filtered_df[filtered_df['OPair Score'] >= score_threshold]

# Filter 3: E-value threshold
filtered_df = filtered_df[filtered_df['Expectation'] <= 0.01]

# Filter 4: Peptides with multiple S/T sites (addresses reviewer concern)
def count_modifiable_sites(peptide_seq):
    """Count S, T, Y residues in peptide"""
    return sum(1 for aa in peptide_seq if aa in ['S', 'T', 'Y'])

filtered_df['n_modifiable_sites'] = filtered_df['Peptide'].apply(count_modifiable_sites)
multi_site_df = filtered_df[filtered_df['n_modifiable_sites'] >= 2]

print(f"PSMs after filtering: {len(multi_site_df)}")
print(f"Unique peptides: {multi_site_df['Peptide'].nunique()}")
```

### Step 1.3: Generate quality score summary table

**Output:** CSV file ranking candidates by combined quality metrics
```python
# Create composite quality score for ranking
multi_site_df['quality_rank'] = (
    multi_site_df['OPair Score'].rank(pct=True) +
    multi_site_df['Hyperscore'].rank(pct=True) +
    (1 - multi_site_df['Expectation'].rank(pct=True))  # lower E-value is better
) / 3

# Sort by quality
ranked_df = multi_site_df.sort_values('quality_rank', ascending=False)

# Save for review
ranked_df.to_csv('ranked_psm_candidates.csv', index=False)
```

---

## Phase 2: Spectral Data Extraction

### Step 2.1: Parse mzML files

**Input:** mzML file(s) containing EThcD spectra

**Packages:**
```python
from pyteomics import mzml
import numpy as np
```

**Code approach:**
```python
def extract_spectrum(mzml_path, scan_number):
    """Extract spectrum by scan number from mzML file"""
    with mzml.read(mzml_path) as reader:
        for spectrum in reader:
            # Check scan number (format varies by instrument)
            spec_id = spectrum['id']
            if f'scan={scan_number}' in spec_id or spectrum.get('index') == scan_number - 1:
                mz_array = spectrum['m/z array']
                intensity_array = spectrum['intensity array']
                precursor_mz = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
                charge = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0].get('charge state', 2)
                return {
                    'mz': mz_array,
                    'intensity': intensity_array,
                    'precursor_mz': precursor_mz,
                    'charge': int(charge),
                    'scan': scan_number
                }
    return None
```

### Step 2.2: Build spectrum index for efficient access

For large mzML files, build an index first:
```python
def build_mzml_index(mzml_path):
    """Build scan number to file offset index"""
    index = {}
    with mzml.read(mzml_path) as reader:
        for spectrum in reader:
            scan_match = re.search(r'scan=(\d+)', spectrum['id'])
            if scan_match:
                scan_num = int(scan_match.group(1))
                index[scan_num] = spectrum
    return index
```

---

## Phase 3: Spectrum Annotation

### Step 3.1: Calculate theoretical fragment ions

**For EThcD spectra, calculate:**
- c-ions (ETD)
- z-ions (ETD)
- b-ions (HCD supplemental activation)
- y-ions (HCD supplemental activation)

**Key consideration:** Include mass shift from O-GlcNAc modification (+203.0794 Da for GlcNAc, but with your workflow it's the azido-HexNAc-PC remnant mass)

```python
from pyteomics import mass, parser

def calculate_theoretical_fragments(peptide_seq, mod_position, mod_mass, charge_states=[1, 2]):
    """
    Calculate theoretical c/z/b/y ions for modified peptide
    
    Parameters:
    - peptide_seq: clean peptide sequence (no modifications)
    - mod_position: 1-based position of modification
    - mod_mass: mass of modification (e.g., 299.1230 for HexNAz-PC)
    - charge_states: list of charge states to consider
    """
    fragments = {'c': [], 'z': [], 'b': [], 'y': []}
    
    n = len(peptide_seq)
    
    for i in range(1, n):
        # N-terminal fragments (b, c ions)
        n_term_seq = peptide_seq[:i]
        n_term_mass = mass.calculate_mass(sequence=n_term_seq)
        
        # Add modification if site is in this fragment
        if mod_position <= i:
            n_term_mass += mod_mass
        
        # b-ion: [M+H]+ - H2O theoretically, but pyteomics handles this
        b_mass = mass.calculate_mass(sequence=n_term_seq, ion_type='b')
        if mod_position <= i:
            b_mass += mod_mass
            
        # c-ion: b + NH3
        c_mass = b_mass + 17.0265
        
        # C-terminal fragments (y, z ions)
        c_term_seq = peptide_seq[i:]
        
        # y-ion
        y_mass = mass.calculate_mass(sequence=c_term_seq, ion_type='y')
        if mod_position > i:
            y_mass += mod_mass
            
        # z-ion: y - NH3 + H (approximately)
        z_mass = y_mass - 16.0187
        
        for z in charge_states:
            fragments['b'].append({'ion': f'b{i}', 'mz': (b_mass + (z-1)*1.007276)/z, 'charge': z})
            fragments['c'].append({'ion': f'c{i}', 'mz': (c_mass + (z-1)*1.007276)/z, 'charge': z})
            fragments['y'].append({'ion': f'y{n-i}', 'mz': (y_mass + (z-1)*1.007276)/z, 'charge': z})
            fragments['z'].append({'ion': f'z{n-i}', 'mz': (z_mass + (z-1)*1.007276)/z, 'charge': z})
    
    return fragments
```

### Step 3.2: Match theoretical to observed peaks

```python
def match_fragments(observed_mz, observed_intensity, theoretical_fragments, tolerance_ppm=20):
    """
    Match observed peaks to theoretical fragments
    
    Returns list of matched ions with their intensities
    """
    matches = []
    
    for ion_type, ions in theoretical_fragments.items():
        for ion in ions:
            theo_mz = ion['mz']
            tol = theo_mz * tolerance_ppm / 1e6
            
            # Find closest observed peak within tolerance
            mz_diff = np.abs(observed_mz - theo_mz)
            min_idx = np.argmin(mz_diff)
            
            if mz_diff[min_idx] <= tol:
                matches.append({
                    'ion_type': ion_type,
                    'ion_label': ion['ion'],
                    'charge': ion['charge'],
                    'theoretical_mz': theo_mz,
                    'observed_mz': observed_mz[min_idx],
                    'intensity': observed_intensity[min_idx],
                    'ppm_error': (observed_mz[min_idx] - theo_mz) / theo_mz * 1e6
                })
    
    return matches
```

---

## Phase 4: Spectrum Visualization

### Step 4.1: Generate annotated spectrum plots

**Recommended package:** `spectrum_utils`

```python
import spectrum_utils.spectrum as sus
import spectrum_utils.plot as sup
import matplotlib.pyplot as plt

def plot_annotated_spectrum(mz, intensity, matches, peptide_seq, mod_site, output_path):
    """
    Generate publication-quality annotated spectrum
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Normalize intensities
    max_int = np.max(intensity)
    norm_intensity = intensity / max_int * 100
    
    # Plot all peaks in gray
    ax.vlines(mz, 0, norm_intensity, colors='gray', linewidth=0.5, alpha=0.5)
    
    # Color scheme for ion types
    colors = {'c': '#1f77b4', 'z': '#ff7f0e', 'b': '#2ca02c', 'y': '#d62728'}
    
    # Plot matched peaks with colors and labels
    for match in matches:
        ion_type = match['ion_type']
        color = colors.get(ion_type, 'black')
        
        ax.vlines(match['observed_mz'], 0, match['intensity']/max_int*100, 
                  colors=color, linewidth=1.5)
        
        # Add ion label
        label = f"{match['ion_label']}"
        if match['charge'] > 1:
            label += f"$^{{{match['charge']}+}}$"
        
        ax.annotate(label, (match['observed_mz'], match['intensity']/max_int*100 + 2),
                    fontsize=8, ha='center', color=color, rotation=90)
    
    # Labels and formatting
    ax.set_xlabel('m/z', fontsize=12)
    ax.set_ylabel('Relative Intensity (%)', fontsize=12)
    ax.set_title(f'{peptide_seq} (O-GlcNAc @ position {mod_site})', fontsize=12)
    
    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], color=colors['c'], label='c ions'),
                       Line2D([0], [0], color=colors['z'], label='z ions'),
                       Line2D([0], [0], color=colors['b'], label='b ions'),
                       Line2D([0], [0], color=colors['y'], label='y ions')]
    ax.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return output_path
```

### Step 4.2: Batch generate spectra for top candidates

```python
def batch_annotate_spectra(ranked_df, mzml_path, output_dir, n_spectra=50):
    """
    Generate annotated spectra for top N candidates
    """
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    results = []
    
    for idx, row in ranked_df.head(n_spectra).iterrows():
        scan_num = row['Scan']  # adjust column name
        peptide = row['Peptide']
        mod_site = row['Modification Site']  # adjust column name
        protein = row['Protein']
        
        # Extract spectrum
        spec_data = extract_spectrum(mzml_path, scan_num)
        if spec_data is None:
            continue
        
        # Calculate theoretical fragments
        # Note: adjust mod_mass based on your actual modification mass
        mod_mass = 299.1230  # HexNAz-PC remnant mass - VERIFY THIS VALUE
        theo_frags = calculate_theoretical_fragments(peptide, mod_site, mod_mass)
        
        # Match fragments
        matches = match_fragments(spec_data['mz'], spec_data['intensity'], theo_frags)
        
        # Calculate coverage metrics
        n_c_ions = len([m for m in matches if m['ion_type'] == 'c'])
        n_z_ions = len([m for m in matches if m['ion_type'] == 'z'])
        coverage = (n_c_ions + n_z_ions) / (2 * (len(peptide) - 1))
        
        # Generate plot
        output_path = os.path.join(output_dir, f'{protein}_{peptide}_{mod_site}_scan{scan_num}.png')
        plot_annotated_spectrum(spec_data['mz'], spec_data['intensity'], 
                                matches, peptide, mod_site, output_path)
        
        results.append({
            'protein': protein,
            'peptide': peptide,
            'mod_site': mod_site,
            'scan': scan_num,
            'n_c_ions': n_c_ions,
            'n_z_ions': n_z_ions,
            'coverage': coverage,
            'n_matches': len(matches),
            'plot_path': output_path
        })
    
    # Save summary
    results_df = pd.DataFrame(results)
    results_df.to_csv(os.path.join(output_dir, 'spectrum_annotation_summary.csv'), index=False)
    
    return results_df
```

---

## Phase 5: Candidate Selection for Figure 6

### Step 5.1: Cross-reference with quantitative data

**Required:** Merge spectral candidates with quantitative log2FC data from your analysis

```python
# Load quantitative data (site-level fold changes)
quant_df = pd.read_csv('site_level_quantification.csv')  # adjust path

# Merge with spectral candidates
merged_df = results_df.merge(quant_df, 
                              left_on=['protein', 'mod_site'],
                              right_on=['Protein', 'Site'],
                              how='inner')
```

### Step 5.2: Add structural annotations

```python
# Load AlphaFold pLDDT scores (from your previous analysis)
plddt_df = pd.read_csv('site_plddt_scores.csv')  # adjust path

# Load UniProt domain annotations
# Option A: Use pre-downloaded domain table
domain_df = pd.read_csv('uniprot_domain_annotations.csv')

# Option B: Query UniProt API
def get_uniprot_domains(accession):
    """Query UniProt for domain annotations"""
    import requests
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    response = requests.get(url)
    if response.ok:
        data = response.json()
        domains = []
        for feature in data.get('features', []):
            if feature['type'] in ['Domain', 'Region', 'Motif', 'Binding site']:
                domains.append({
                    'type': feature['type'],
                    'description': feature.get('description', ''),
                    'start': feature['location']['start']['value'],
                    'end': feature['location']['end']['value']
                })
        return domains
    return []

# Check if modification site falls within a domain
def site_in_domain(site_position, domains):
    """Check if site falls within any annotated domain"""
    for domain in domains:
        if domain['start'] <= site_position <= domain['end']:
            return domain['description']
    return None
```

### Step 5.3: Filter for Figure 6 candidates

**For structured region example (like current HYOU1):**
```python
structured_candidates = merged_df[
    (merged_df['pLDDT'] >= 70) &  # high confidence = structured
    (merged_df['log2FC'] > 1) &    # upregulated
    (merged_df['coverage'] >= 0.3) &  # good spectral coverage
    (merged_df['domain'].notna())  # has annotated domain
].sort_values('coverage', ascending=False)
```

**For IDR example (like current HOXA13):**
```python
idr_candidates = merged_df[
    (merged_df['pLDDT'] < 50) &    # low confidence = IDR
    (merged_df['log2FC'].abs() < 0.5) &  # minimal change
    (merged_df['coverage'] >= 0.3)
].sort_values('coverage', ascending=False)
```

---

## Output Files

| File | Description |
|------|-------------|
| `ranked_psm_candidates.csv` | All PSMs ranked by quality score |
| `annotated_spectra/` | Directory of annotated spectrum PNGs |
| `spectrum_annotation_summary.csv` | Summary metrics for each annotated spectrum |
| `figure6_candidates.csv` | Final filtered candidates for Figure 6 |

---

## User Input Required

Before execution, please provide:

1. **File paths:**
   - [ ] O-Pair PSM output file path
   - [ ] mzML file path(s)
   - [ ] Site-level quantification table (with log2FC values)
   - [ ] pLDDT annotation table (from AlphaFold analysis)

2. **Modification mass:**
   - [ ] Confirm exact mass of HexNAz-PC remnant in your search (299.1230 Da?)
   - [ ] Confirm mass of GAO-modified form if applicable

3. **Column name mappings:**
   - [ ] Scan number column name in PSM file
   - [ ] Peptide sequence column name
   - [ ] Modification site column name
   - [ ] Protein accession column name

4. **Filtering thresholds:**
   - [ ] Minimum O-Pair score (or use percentile-based)
   - [ ] pLDDT cutoff for structured (≥70?) vs IDR (<50?)
   - [ ] log2FC threshold for "upregulated" vs "unchanged"

---

## Dependencies

```bash
pip install pyteomics spectrum_utils pandas numpy matplotlib requests
```

---

## Notes

- The modification mass values in this plan are placeholders - verify against your actual search parameters
- Fragment ion calculations assume standard amino acid masses - check if your workflow includes TMT labels
- For TMT-labeled peptides, add TMT mass (+229.1629) to N-terminus and lysines
- Consider adding neutral loss peaks (e.g., HexNAc loss) if reviewing HCD spectra
