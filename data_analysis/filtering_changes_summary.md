# Summary: O-GlcNAc Site Filtering Changes for Manuscript Revision

## Overview

This document summarizes the filtering changes applied to O-GlcNAc site data for manuscript revision. The key change is applying a **probability threshold of 0.75** to Level1b PSMs to ensure confident site localization.

---

## Filtering Criteria

### Previous Approach
All Level1 and Level1b PSMs were included regardless of site probability.

### Updated Approach
Apply a probability threshold of 0.75 to Level1b PSMs:

| Confidence Level | Description | Filtering Rule |
|------------------|-------------|----------------|
| **Level1** | EThcD spectra with c/z ion evidence | Keep all |
| **Level1b** | HCD-only spectra | Keep only if probability ≥ 0.75 |

---

## Rationale

1. **Level1 PSMs (EThcD)**: All 1,872 Level1 PSMs already have probability ≥ 0.75 (range: 0.751–1.000, mean: 0.931). The O-Pair algorithm inherently assigns high confidence to spectra with good c/z ion coverage.

2. **Level1b PSMs (HCD-only)**: These lack c/z ions for site localization. PSMs with probability < 0.75 have weaker localization evidence. A total of 471 Level1b PSMs have probability < 0.75 and are removed.

3. **No contradictions**: Manual analysis of b/y ions in the removed 471 PSMs showed that 100% of assigned sites are consistent with the spectral evidence—none were contradicted.

---

## Impact on PSM Counts (Table S6)

| Cell Type | Before Filtering | After Filtering | Removed (Level1b prob < 0.75) |
|-----------|------------------|-----------------|-------------------------------|
| HEK293T | 1,045 | 832 | 213 |
| HepG2 | 599 | 508 | 91 |
| Jurkat | 730 | 563 | 167 |
| **Total** | **2,374** | **1,903** | **471 (19.8%)** |

---

## Impact on Unique Site Counts (Table S7)

| Cell Type | Sites Before Filtering | Sites After Filtering | Removed Sites |
|-----------|------------------------|----------------------|---------------|
| HEK293T | 366 | 313 | 53 |
| HepG2 | 290 | 257 | 33 |
| Jurkat | 283 | 226 | 57 |
| **Total** | **939** | **796** | **143 (15.2%)** |

Sites used in regression analysis (with complete structural features):
- HEK293T: 277
- HepG2: 217
- Jurkat: 195
- **Total: 689 sites**

---

## Validation: b/y Ion Localization Analysis

We performed a detailed analysis of all 471 removed Level1b PSMs (prob < 0.75) using HCD b/y ion evidence:

| Category | Count | % of Removed |
|----------|-------|--------------|
| **Localized** (uniquely determined by b/y ions) | 254 | 53.9% |
| **Not Localized** (multiple sites possible) | 217 | 46.1% |
| **Contradicted** (assigned site ruled out) | 0 | 0% |

### Annotation Statistics
- Average sequence coverage: 58.5%
- False match rate (FMR): 0.4%
- Site-determining b ions: 1,388 (66.3%)
- Site-determining y ions: 707 (33.7%)

### Key Finding
In **100% of cases**, the assigned glycosylation site is consistent with b/y ion evidence. The removed PSMs are **not necessarily wrong**—they simply lack the c/z ion evidence required by the O-Pair algorithm for high confidence scoring.

---

## Validation: Comparison of Filtering Approaches

We compared two filtering approaches:
- **Option 1 (Site-level filtering)**: Quantify all PSMs, then keep sites with at least one high-confidence PSM
- **Option 2 (PSM-level filtering)**: Filter PSMs first, then quantify only high-confidence PSMs

### Results

| Metric | Value |
|--------|-------|
| Total sites compared | 911 |
| Overall logFC correlation | **0.9875** |
| Mean |logFC difference| | **0.053** |
| Sites with |difference| > 0.5 | 10 (1.1%) |
| Sites with |difference| > 1.0 | 7 (0.8%) |

### Significance Changes (adj.P.Val < 0.05)

| Category | Count | % |
|----------|-------|---|
| Both significant | 478 | 52.5% |
| Both non-significant | 412 | 45.2% |
| Gained significance | 10 | 1.1% |
| Lost significance | 11 | 1.2% |

### Conclusion
The two approaches produce nearly identical results (r = 0.9875). Site-level filtering (Option 1) is used because it retains valid intensity data from all PSMs while ensuring sites have high-confidence localization evidence.

---

## Effect on Regression Analysis (Figure 6)

| Metric | Before Filtering | After Filtering |
|--------|------------------|-----------------|
| Total sites (n) | 939 | 689 |
| R-squared | ~0.10 | 0.099 |
| Adj R-squared | ~0.08 | 0.082 |

### Significant Predictors (Unchanged)
| Predictor | Standardized β | p-value |
|-----------|----------------|---------|
| IDR (vs Structured) | -0.82 | 3.1×10⁻⁷ *** |
| pPSE_24 (Full Sphere) | -0.20 | 0.0025 ** |

The filtering reduces sample size by ~27% but does not change the regression conclusions.

---

## Updated Supporting Tables

### Table S6: Identification of O-GlcNAcylation Sites
- **Content**: PSM-level identification data
- **Filtering applied**: Removed Level1b PSMs with probability < 0.75
- **Final counts**: 1,903 PSMs (832 HEK293T, 508 HepG2, 563 Jurkat)

### Table S7: Abundance Changes of O-GlcNAcylation Sites
- **Content**: Site-level differential expression (logFC, adj.P.Val)
- **Filtering applied**: Removed sites with no high-confidence PSMs
- **Final counts**: 796 sites (313 HEK293T, 257 HepG2, 226 Jurkat)

---

## Updated Figure 6

Figure 6 panels A-D have been regenerated using the filtered data:
- **Figure6A**: logFC distribution by cell type
- **Figure6B**: Secondary structure distribution (donut plot)
- **Figure6C**: Standardized regression coefficients
- **Figure6D**: IDR vs Structured effect by cell type

Output directory: `Figures/Figure6_filtered/`

---

## Suggested Manuscript Text

### Methods Section

> **Site localization filtering.** For site-level analysis, we applied a probability threshold of 0.75 to Level1b PSMs (HCD-only spectra) to ensure confident site localization. Level1 PSMs (EThcD spectra with c/z ion evidence) were retained without additional filtering, as all met the probability threshold (range: 0.751–1.000). This filtering removed 471 of 2,374 PSMs (19.8%), resulting in 1,903 high-confidence PSMs representing 796 unique O-GlcNAc sites for downstream quantitative analysis.

### Results Section (if discussing filtering)

> To ensure confident site localization, we applied a probability threshold of 0.75 to Level1b identifications, removing 471 PSMs (19.8%) with lower localization confidence. Analysis of the removed spectra using b/y ion evidence confirmed that none of the assigned sites were contradicted by the spectral data, with 53.9% uniquely localizable by b/y ions alone. The remaining 46.1% had multiple possible sites but the assigned position remained consistent with all observed ions.

### Supplementary Information

> **Site probability filtering.** O-GlcNAc site identifications were filtered to retain only high-confidence localizations. Level1 identifications (EThcD spectra with c/z ion evidence) were retained in full, as all had site probabilities ≥ 0.75 (range: 0.751–1.000). Level1b identifications (HCD-only spectra) were filtered to retain only those with site probability ≥ 0.75. This conservative approach removed 471 Level1b PSMs (19.8% of total), resulting in 1,903 PSMs across three cell types. Detailed analysis of the removed PSMs using b/y ion constraints demonstrated that while these sites lacked sufficient evidence for high O-Pair probability scores, the assigned localizations were not contradicted by the available spectral evidence.

---

## Files Generated

| File | Description |
|------|-------------|
| `supporting_table_S6.xlsx` | O-GlcNAc site identification (filtered) |
| `supporting_table_S7.xlsx` | O-GlcNAc site abundance changes (filtered) |
| `OGlcNAc_site_features_filtered.csv` | Site features for regression analysis |
| `Figure6_filtered/*.pdf` | Updated Figure 6 panels |
| `filtering_comparison/` | Comparison of filtering approaches |
| `low_prob_spectra/` | Analysis of removed Level1b PSMs |

---

## Summary Statistics for Quick Reference

| Metric | Value |
|--------|-------|
| Probability threshold | 0.75 |
| PSMs removed | 471 (19.8%) |
| PSMs retained | 1,903 |
| Unique sites removed | 143 (15.2%) |
| Unique sites retained | 796 |
| Sites for regression | 689 |
| logFC correlation (Option 1 vs 2) | 0.9875 |
| Removed PSMs with b/y localization | 53.9% |
| Removed PSMs contradicted | 0% |
