# Figure 6 Filtered Analysis Summary

## Filtering Criteria

Removed Level1b PSMs with site probability < 0.75 (low confidence HCD-based localizations).

- **Keep**: Level1 (all)
- **Keep**: Level1b with probability >= 0.75
- **Remove**: Level1b with probability < 0.75

## Filtering Results

### PSM Counts

| Cell Type | Before Filtering | After Filtering | Removed |
|-----------|------------------|-----------------|---------|
| HEK293T | 1,045 | 832 | 213 |
| HepG2 | 599 | 508 | 91 |
| Jurkat | 730 | 563 | 167 |
| **Total** | **2,374** | **1,903** | **471** |

### Unique Sites for Differential Analysis

| Cell Type | Sites After Filtering |
|-----------|----------------------|
| HEK293T | 313 |
| HepG2 | 257 |
| Jurkat | 226 |
| **Total (per cell type sum)** | **796** |
| **Total (unique across all)** | **559** |

Note: 796 is the sum of sites per cell type (used for regression analysis where each cell type is analyzed separately). 559 is the number of unique site identities across all cell types, accounting for overlapping sites.

### Site Overlaps Across Cell Types

| Overlap | Count |
|---------|-------|
| HEK293T ∩ HepG2 | 113 |
| HEK293T ∩ Jurkat | 107 |
| HepG2 ∩ Jurkat | 94 |
| All three cell types | 77 |

## Regression Analysis Results (Filtered Data)

### Model Summary

- **Sample size**: n = 689 sites with complete data
  - HEK293T: 277
  - HepG2: 217
  - Jurkat: 195
- **R-squared**: 0.0989
- **Adjusted R-squared**: 0.0816

### Standardized Coefficients

| Predictor | Std. Estimate | p-value | Significance |
|-----------|---------------|---------|--------------|
| IDR (vs Structured) | -0.823 | 3.09e-07 | *** |
| pPSE_24 (Full Sphere) | -0.195 | 0.0025 | ** |
| Serine (vs Thr) | 0.094 | 0.112 | |
| Jurkat (vs HEK293T) | -0.066 | 0.328 | |
| pPSE_12 (Side Chain) | 0.064 | 0.137 | |
| pLDDT (Confidence) | 0.060 | 0.219 | |
| pI (7-mer) | -0.054 | 0.058 | . |
| Hydrophobicity | 0.036 | 0.196 | |
| Sites per Protein | -0.027 | 0.371 | |
| HepG2 (vs HEK293T) | 0.024 | 0.716 | |
| Protein Abundance | 0.015 | 0.617 | |

Significance codes: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.1

## Key Finding

**The trend is preserved after filtering.**

The IDR effect remains highly significant (p = 3.09e-07):

- Sites in **structured regions** show **larger positive logFC** (increased O-GlcNAcylation with tunicamycin treatment)
- Sites in **IDR regions** show **smaller/negative logFC** (less change or decreased O-GlcNAcylation)

This confirms that the structural context (IDR vs structured) association with tunicamycin response is robust and not driven by low-confidence site localizations.

\newpage

## Verification of Figure 6E/6F Example Sites

The final examples for Figure 6E (structured region) and Figure 6F (IDR region) were verified to ensure they are not affected by the filtering.

### PRDX6 Y89 (Figure 6E - Structured Region)

| Cell Type | Confidence | Probability | Status |
|-----------|------------|-------------|--------|
| HEK293T | Level1b | 0.500 | Removed |
| HEK293T | **Level1** | **0.936** | **Kept (EThcD)** |
| Jurkat | **Level1** | **0.776** | **Kept (EThcD)** |
| Jurkat | **Level1** | **0.983** | **Kept (EThcD)** |
| Jurkat | Level1b | 0.500 | Removed |

**Differential Expression (HEK293T):**

- log2FC: **1.5002** (FC = 2.83)
- adj.P.Val: **0.00054** (significant)
- Status: Upregulated in structured region

### EWSR1 S274 (Figure 6F - IDR Region)

| Cell Type | Confidence | Probability | Status |
|-----------|------------|-------------|--------|
| HEK293T | **Level1** | **0.983** | **Kept (EThcD)** |

**Differential Expression (HEK293T):**

- log2FC: **-0.1492** (FC = 0.90)
- adj.P.Val: 0.507 (not significant)
- Status: Essentially unchanged in IDR region

### Verification Summary

| Site | Expected FC | Expected log2FC | Observed log2FC | Match |
|------|-------------|-----------------|-----------------|-------|
| PRDX6 Y89 | 2.83 | 1.5008 | 1.5002 | Yes |
| EWSR1 S274 | 0.90 | -0.1520 | -0.1492 | Yes |

**Conclusion:** Both examples have Level1 (EThcD) evidence and are not affected by the filtering. The fold changes match expected values, confirming the examples properly illustrate the key finding.

\newpage

## Figure 6A: logFC Distribution by Cell Type

![Figure 6A (Filtered): logFC distribution showing the fold change distribution across HEK293T, HepG2, and Jurkat cells after removing low probability Level1b PSMs.](Figure6A_filtered.png){width=80%}

\newpage

## Figure 6B: Secondary Structure Distribution

![Figure 6B (Filtered): Donut plot showing the distribution of O-GlcNAc sites across secondary structure categories (Strand, Bend, Unstructured).](Figure6B_filtered.png){width=60%}

\newpage

## Figure 6C: Standardized Regression Coefficients

![Figure 6C (Filtered): Bar plot of standardized regression coefficients showing the effect of each predictor (1 SD change) on log2FC. Significant predictors (p < 0.05) are highlighted.](Figure6C_filtered.png){width=90%}

\newpage

## Figure 6D: IDR vs Structured Effect Across Cell Types

![Figure 6D (Filtered): Mean log2FC comparison between IDR and Structured regions across the three cell types. The consistent pattern shows structured regions have higher mean logFC than IDR regions.](Figure6D_filtered.png){width=60%}
