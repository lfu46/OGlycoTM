# O-GlcNAc Site Structural Regression Analysis Plan

## Project Overview

### Objective
Identify which structural and sequence features are associated with treatment-induced 
fold changes at O-GlcNAc modification sites using multivariate linear regression.

### Scientific Question
When cells are treated with a perturbation (e.g., OGT/OGA inhibitor, metabolic stress, 
signaling pathway activation), O-GlcNAc levels change at specific sites. We want to 
understand: **What structural or sequence characteristics make a site more or less 
responsive to treatment?**

### Background and Rationale
This analysis approach is informed by Potel et al. 2025 (Nature Structural & Molecular 
Biology, "Deep quantitative glycoprofiling"), which used structural features from 
AlphaFold predictions to understand glycosylation site properties. Key insights from 
that work:

1. **pPSE (prediction-aware part-sphere exposure)** is superior to binary "buried vs 
   exposed" classification because it provides continuous measurement of solvent 
   accessibility while accounting for AlphaFold prediction confidence.

2. **IDR (intrinsically disordered region)** status captures information distinct from 
   accessibility—a site can be exposed but in a disordered region, or buried within 
   an ordered domain.

3. **Multiple sites per protein are not independent**—they share the same cellular 
   environment, protein turnover rate, and trafficking. Statistical analysis must 
   account for this clustering.

4. **O-GlcNAc has unique biology** compared to N-glycosylation: it occurs on Ser/Thr 
   residues in cytoplasm/nucleus, is catalyzed by a single enzyme pair (OGT/OGA), 
   and exhibits extensive crosstalk with phosphorylation.

---

## Feature Specification

### 1. Structural Accessibility Features

These features quantify how exposed the modification site is to the solvent and 
potentially to modifying enzymes (OGT/OGA).

| Feature | Type | Source | Description |
|---------|------|--------|-------------|
| `pPSE_24` | Continuous | structuremap + AlphaFold | Part-sphere exposure with 24 Å radius and 180° angle. Counts residues in structural proximity, weighted by AlphaFold confidence. Higher values = more buried/ordered. |
| `pPSE_ratio` | Continuous | Derived | Ratio of pPSE_12 to pPSE_24. Captures local vs global exposure. Calculate as: `pPSE_12 / pPSE_24`. |

**Note on pPSE calculation:**
- Use the `structuremap` Python package (version 0.0.9 or later)
- Requires AlphaFold structure files (.pdb or .cif format)
- pPSE_12 uses 12 Å radius and 70° angle (local environment)
- pPSE_24 uses 24 Å radius and 180° angle (global environment)

**Why use ratio instead of both raw values:**
pPSE_12 and pPSE_24 are often highly correlated (>0.7). Including both causes 
multicollinearity. The ratio captures whether a site has different local vs global 
exposure characteristics.

### 2. Secondary Structure Features

| Feature | Type | Source | Description |
|---------|------|--------|-------------|
| `secondary_structure` | Categorical | structuremap + AlphaFold | Categories: helix, strand, turn, bend, coil |

**For regression, convert to dummy variables:**
- `is_helix`: 1 if α-helix, 0 otherwise
- `is_strand`: 1 if β-strand, 0 otherwise  
- `is_turn`: 1 if turn, 0 otherwise
- `is_bend`: 1 if bend, 0 otherwise
- Reference category: coil (most common, will be absorbed into intercept)

### 3. Disorder Feature

| Feature | Type | Source | Description |
|---------|------|--------|-------------|
| `is_IDR` | Binary | Derived from pPSE or pLDDT | Whether site is in intrinsically disordered region |

**Definition options (use one consistently):**
- Option A: `pPSE_24_smoothed < 34.27` (threshold from structuremap paper)
- Option B: `pLDDT < 50` (low AlphaFold confidence indicates disorder)
- Option C: Use MobiDB or DisProt database annotations

### 4. Domain Context Feature

| Feature | Type | Source | Description |
|---------|------|--------|-------------|
| `in_domain` | Binary | InterPro/UniProt | Whether site falls within an annotated protein domain |

**Data source:**
- Download domain annotations from InterPro (https://www.ebi.ac.uk/interpro/)
- Or use UniProt feature annotations (DOMAIN, REGION)
- Check if site position falls within any domain boundaries

### 5. Sequence Property Features

Calculate these from the 7-mer peptide sequence centered on the modification site 
(3 residues on each side of the modified Ser/Thr).

| Feature | Type | Calculation | Description |
|---------|------|-------------|-------------|
| `pI_7mer` | Continuous | BioPython ProteinAnalysis | Isoelectric point of flanking sequence. Reflects local charge environment. |
| `hydrophobicity_7mer` | Continuous | GRAVY score | Grand average of hydropathy. Positive = hydrophobic, negative = hydrophilic. |
| `is_serine` | Binary | Direct from data | 1 if modified residue is Serine, 0 if Threonine. OGT may have different kinetics for Ser vs Thr. |

**Calculation code pattern:**
```python
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def calculate_sequence_features(sequence_7mer):
    """Calculate pI and hydrophobicity for a 7-mer sequence."""
    # Handle non-standard amino acids if present
    clean_seq = ''.join([aa for aa in sequence_7mer if aa in 'ACDEFGHIKLMNPQRSTVWY'])
    
    if len(clean_seq) < 3:
        return np.nan, np.nan
    
    analysis = ProteinAnalysis(clean_seq)
    pI = analysis.isoelectric_point()
    hydrophobicity = analysis.gravy()
    
    return pI, hydrophobicity
```

### 6. PTM Crosstalk Feature

| Feature | Type | Source | Description |
|---------|------|--------|-------------|
| `near_phosphosite` | Binary | PhosphoSitePlus | Whether a known phosphorylation site exists within ±5 residues |

**Rationale:**
O-GlcNAc and phosphorylation exhibit extensive crosstalk—they can compete for the 
same or adjacent Ser/Thr residues, and their regulatory enzymes can influence each 
other. Sites near phosphosites may respond differently to perturbation.

**Data source:**
- Download from PhosphoSitePlus (https://www.phosphosite.org/)
- Filter for organism and high-confidence sites
- Check if any phosphosite falls within window [site_position - 5, site_position + 5]

### 7. Protein-Level Covariates

These are not features of biological interest but confounders that must be controlled.

| Feature | Type | Source | Description |
|---------|------|--------|-------------|
| `log_protein_abundance` | Continuous | Total proteome data | Log2 protein intensity from unenriched proteome. More abundant proteins have more reliable quantification. |
| `sites_per_protein` | Continuous | Derived | Number of O-GlcNAc sites detected on this protein. Proteins with many sites may show different regulation patterns (competition for OGT). |

---

## Data Requirements

### Input Data Files

1. **Glycopeptide quantification file** (e.g., from MSFragger-Glyco, Byonic, or similar)
   - Required columns:
     - `protein_id`: UniProt accession
     - `site_position`: Position of modified residue in protein sequence
     - `modified_residue`: S or T
     - `peptide_sequence`: Full peptide sequence
     - `fold_change`: Log2 fold change (treatment / control)
     - `adj_pvalue`: Adjusted p-value from differential analysis
   
2. **Protein abundance file** (from total proteome analysis)
   - Required columns:
     - `protein_id`: UniProt accession
     - `abundance`: Protein intensity or normalized abundance

3. **AlphaFold structure files**
   - One .pdb or .cif file per protein
   - Download from AlphaFold Database (https://alphafold.ebi.ac.uk/)
   - Naming convention: AF-{UniProt_ID}-F1-model_v4.pdb

4. **UniProt annotations**
   - Protein sequences (for extracting 7-mer)
   - Domain boundaries

5. **PhosphoSitePlus phosphorylation sites**
   - Filtered for organism of interest

### Output Data Structure

After feature extraction, the analysis-ready dataframe should have this structure:

| Column | Type | Example |
|--------|------|---------|
| protein_id | string | "P12345" |
| site_position | int | 245 |
| peptide_sequence | string | "PEPTIDES*EQUENCE" |
| fold_change | float | 1.35 |
| adj_pvalue | float | 0.023 |
| pPSE_24 | float | 156.3 |
| pPSE_ratio | float | 0.42 |
| is_helix | int | 0 |
| is_strand | int | 1 |
| is_turn | int | 0 |
| is_bend | int | 0 |
| is_IDR | int | 0 |
| in_domain | int | 1 |
| pI_7mer | float | 5.82 |
| hydrophobicity_7mer | float | -0.34 |
| is_serine | int | 1 |
| near_phosphosite | int | 1 |
| log_protein_abundance | float | 25.6 |
| sites_per_protein | int | 3 |

---

## Statistical Analysis Approach

### Primary Model: OLS with Clustered Standard Errors

**Model formula:**
```
fold_change ~ pPSE_24 + pPSE_ratio + 
              is_helix + is_strand + is_turn + is_bend +
              is_IDR + in_domain + 
              pI_7mer + hydrophobicity_7mer + is_serine +
              near_phosphosite + 
              log_protein_abundance + sites_per_protein
```

**Why clustered standard errors:**
Multiple O-GlcNAc sites on the same protein are not independent observations. They 
share protein abundance, turnover rate, subcellular localization, and cellular context. 
Standard OLS assumes independence, which would underestimate standard errors and 
inflate significance. Clustering by `protein_id` corrects for this.

**Implementation:**
```python
import statsmodels.formula.api as smf

model = smf.ols(formula, data=df).fit(
    cov_type='cluster',
    cov_kwds={'groups': df['protein_id']}
)
```

### Alternative Model: Mixed Effects

If clustered standard errors produce unstable estimates (e.g., many proteins with 
only 1 site), use a mixed effects model with protein as random intercept:
```python
model_mixed = smf.mixedlm(
    formula,
    data=df,
    groups=df['protein_id']
).fit()
```

### Variable Standardization

Standardize all continuous variables before fitting to:
1. Make coefficients comparable (all on same scale)
2. Improve numerical stability
3. Allow interpretation as "effect of 1 SD change"
```python
from sklearn.preprocessing import StandardScaler

continuous_vars = ['pPSE_24', 'pPSE_ratio', 'pI_7mer', 
                   'hydrophobicity_7mer', 'log_protein_abundance', 'sites_per_protein']

scaler = StandardScaler()
df[continuous_vars] = scaler.fit_transform(df[continuous_vars])

# Save scaler parameters for reporting
scaling_params = pd.DataFrame({
    'variable': continuous_vars,
    'mean': scaler.mean_,
    'std': scaler.scale_
})
```

### Multicollinearity Check

Before interpreting coefficients, verify no severe multicollinearity:
```python
from statsmodels.stats.outliers_influence import variance_inflation_factor

def calculate_vif(df, features):
    """Calculate Variance Inflation Factor for each feature."""
    X = df[features].copy()
    X = sm.add_constant(X)
    
    vif_data = pd.DataFrame()
    vif_data['feature'] = X.columns
    vif_data['VIF'] = [
        variance_inflation_factor(X.values, i) 
        for i in range(X.shape[1])
    ]
    return vif_data

vif_results = calculate_vif(df, feature_columns)
print(vif_results)

# Flag features with VIF > 5 (some use threshold of 10)
problematic = vif_results[vif_results['VIF'] > 5]
if len(problematic) > 0:
    print("WARNING: High multicollinearity detected:")
    print(problematic)
```

**If multicollinearity detected:**
- Check correlation matrix to identify which features are correlated
- Consider dropping one of the correlated features
- Or combine them (e.g., PCA, ratio)

### Residual Diagnostics

Check model assumptions by examining residuals:
```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# 1. Residuals vs Fitted (check for patterns)
axes[0,0].scatter(model.fittedvalues, model.resid, alpha=0.3)
axes[0,0].axhline(y=0, color='red', linestyle='--')
axes[0,0].set_xlabel('Fitted values')
axes[0,0].set_ylabel('Residuals')
axes[0,0].set_title('Residuals vs Fitted')

# 2. Q-Q plot (check normality of residuals)
sm.qqplot(model.resid, line='45', ax=axes[0,1])
axes[0,1].set_title('Q-Q Plot')

# 3. Histogram of residuals
axes[1,0].hist(model.resid, bins=50, edgecolor='black')
axes[1,0].set_xlabel('Residuals')
axes[1,0].set_title('Distribution of Residuals')

# 4. Residuals vs key predictor (check linearity assumption)
axes[1,1].scatter(df['pPSE_24'], model.resid, alpha=0.3)
axes[1,1].axhline(y=0, color='red', linestyle='--')
axes[1,1].set_xlabel('pPSE_24')
axes[1,1].set_ylabel('Residuals')
axes[1,1].set_title('Residuals vs pPSE_24')

plt.tight_layout()
plt.savefig('residual_diagnostics.png', dpi=150)
```

---

## Expected Outputs

### 1. Coefficient Table (CSV)

| feature | coefficient | std_error | t_statistic | p_value | ci_lower | ci_upper | significant |
|---------|-------------|-----------|-------------|---------|----------|----------|-------------|
| pPSE_24 | 0.152 | 0.034 | 4.47 | 0.00001 | 0.085 | 0.219 | TRUE |
| ... | ... | ... | ... | ... | ... | ... | ... |

### 2. Forest Plot (PNG)

Horizontal bar plot showing:
- Coefficient point estimate for each feature
- 95% confidence interval as error bars
- Vertical line at zero (no effect)
- Features sorted by absolute coefficient magnitude
- Color coding: significant (p < 0.05) vs non-significant
```python
def create_forest_plot(model, output_path='forest_plot.png'):
    """Create forest plot of regression coefficients."""
    # Extract results (exclude intercept)
    results = pd.DataFrame({
        'feature': model.params.index[1:],
        'coef': model.params.values[1:],
        'ci_lower': model.conf_int()[0].values[1:],
        'ci_upper': model.conf_int()[1].values[1:],
        'pvalue': model.pvalues.values[1:]
    })
    
    # Sort by absolute coefficient
    results['abs_coef'] = results['coef'].abs()
    results = results.sort_values('abs_coef', ascending=True)
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 8))
    
    y_pos = range(len(results))
    colors = ['#d62728' if p < 0.05 else '#7f7f7f' for p in results['pvalue']]
    
    ax.barh(y_pos, results['coef'], color=colors, alpha=0.7)
    ax.errorbar(results['coef'], y_pos,
                xerr=[results['coef'] - results['ci_lower'],
                      results['ci_upper'] - results['coef']],
                fmt='none', color='black', capsize=3)
    
    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(results['feature'])
    ax.set_xlabel('Standardized Coefficient (95% CI)')
    ax.set_title('Feature Associations with Fold Change')
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#d62728', alpha=0.7, label='p < 0.05'),
        Patch(facecolor='#7f7f7f', alpha=0.7, label='p ≥ 0.05')
    ]
    ax.legend(handles=legend_elements, loc='lower right')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    
    return results
```

### 3. VIF Table (CSV)

| feature | VIF | status |
|---------|-----|--------|
| pPSE_24 | 1.82 | OK |
| pPSE_ratio | 2.15 | OK |
| ... | ... | ... |

### 4. Model Summary (TXT)

Full statsmodels summary output including:
- R-squared and adjusted R-squared
- F-statistic and p-value
- Number of observations
- All coefficient statistics

### 5. Scaling Parameters (CSV)

For reproducing standardization and back-transforming if needed:

| variable | mean | std |
|----------|------|-----|
| pPSE_24 | 145.3 | 52.7 |
| ... | ... | ... |

---

## Interpretation Guidelines

### Coefficient Interpretation

After standardization, coefficients represent:

> "A 1 standard deviation increase in [feature] is associated with a [coefficient] 
> change in log2 fold change, holding all other features constant."

**Example:**
If coefficient for `pPSE_24` is 0.15:
- Sites with higher pPSE_24 (more buried/structured) show larger fold changes
- Specifically: sites 1 SD above mean pPSE_24 have fold changes 0.15 log2 units 
  higher than sites at mean pPSE_24

### For Binary Features

Coefficients represent:
> "The difference in log2 fold change between sites with feature=1 versus feature=0, 
> holding all other features constant."

**Example:**
If coefficient for `near_phosphosite` is -0.25:
- Sites near phosphosites show fold changes 0.25 log2 units lower than sites not 
  near phosphosites
- This could indicate competition or regulatory crosstalk

### Comparing Feature Importance

Because continuous variables are standardized, coefficient magnitudes are directly 
comparable. Larger absolute coefficients indicate stronger associations.

For comparing continuous vs binary features, note that binary features represent 
discrete group differences, while continuous features represent per-SD effects.

### Caveats

1. **Association, not causation**: Regression identifies correlations; it cannot prove 
   that a structural feature causes differential response.

2. **Unmeasured confounders**: Other factors (e.g., local OGT concentration, chromatin 
   accessibility for nuclear sites) may explain associations.

3. **Model assumptions**: Results are valid only if linearity, homoscedasticity, and 
   normality assumptions are reasonably met. Check residual diagnostics.

4. **Multiple testing**: With many features, some may appear significant by chance. 
   Consider the overall pattern rather than individual p-values.

---

## Code Organization

### Suggested File Structure
```
project/
├── docs/
│   └── analysis_plan.md              # This document
├── data/
│   ├── raw/
│   │   ├── glycopeptides.csv         # Input quantification data
│   │   ├── proteome.csv              # Protein abundance data
│   │   └── phosphosites.csv          # PhosphoSitePlus data
│   ├── external/
│   │   ├── alphafold_structures/     # Downloaded AlphaFold PDBs
│   │   └── uniprot_annotations.tsv   # Domain annotations
│   └── processed/
│       └── features_matrix.csv       # Analysis-ready data
├── src/
│   ├── feature_extraction.py         # Functions to calculate features
│   ├── regression_analysis.py        # Main analysis script
│   └── visualization.py              # Plotting functions
├── results/
│   ├── tables/
│   │   ├── coefficients.csv
│   │   ├── vif_table.csv
│   │   └── scaling_params.csv
│   ├── figures/
│   │   ├── forest_plot.png
│   │   └── residual_diagnostics.png
│   └── model_summary.txt
└── requirements.txt
```

### Key Dependencies
```
# requirements.txt
pandas>=1.5.0
numpy>=1.23.0
statsmodels>=0.13.0
scikit-learn>=1.1.0
biopython>=1.79
matplotlib>=3.6.0
seaborn>=0.12.0
structuremap>=0.0.9
```

---

## Implementation Checklist

### Phase 1: Data Preparation
- [ ] Load glycopeptide quantification data
- [ ] Load protein abundance data
- [ ] Download AlphaFold structures for all proteins in dataset
- [ ] Download UniProt annotations (domains)
- [ ] Download PhosphoSitePlus phosphorylation sites

### Phase 2: Feature Extraction
- [ ] Calculate pPSE_24 and pPSE_12 using structuremap
- [ ] Calculate pPSE_ratio
- [ ] Extract secondary structure annotations
- [ ] Determine IDR status for each site
- [ ] Map domain annotations to sites
- [ ] Extract 7-mer sequences and calculate pI and hydrophobicity
- [ ] Identify modified residue type (Ser/Thr)
- [ ] Map phosphosites and calculate proximity
- [ ] Calculate sites_per_protein
- [ ] Merge all features into single dataframe

### Phase 3: Data Quality
- [ ] Check for missing values in features
- [ ] Examine feature distributions
- [ ] Calculate correlation matrix
- [ ] Flag potential outliers

### Phase 4: Statistical Analysis
- [ ] Standardize continuous variables
- [ ] Calculate VIF for multicollinearity check
- [ ] Fit OLS model with clustered standard errors
- [ ] Generate residual diagnostic plots
- [ ] Extract coefficient table with confidence intervals

### Phase 5: Output Generation
- [ ] Save coefficient table as CSV
- [ ] Create forest plot
- [ ] Save VIF table
- [ ] Save model summary
- [ ] Save scaling parameters

---

## Contact and Questions

If any aspect of this analysis plan is unclear or needs modification, please ask for 
clarification before implementing. Key decision points that may need discussion:

1. Choice of IDR definition method
2. Handling of proteins without AlphaFold structures
3. Treatment of sites with missing feature values
4. Threshold for VIF (5 vs 10)
5. Whether to use clustered SE vs mixed effects model