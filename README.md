# Systematic Quantification of Protein *O*-GlcNAcylation under *N*-Glycosylation Inhibition

[![Published in Analytical Chemistry](https://img.shields.io/badge/Anal.%20Chem.-2026%2C%2098%2C%2015689–15699-B31B1B)](https://doi.org/10.1021/acs.analchem.6c00972)
[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.analchem.6c00972-blue)](https://doi.org/10.1021/acs.analchem.6c00972)
[![Data: PXD073249](https://img.shields.io/badge/ProteomeXchange-PXD073249-1f9c5a)](https://www.ebi.ac.uk/pride/archive/projects/PXD073249)
[![License: MIT](https://img.shields.io/badge/Code%20License-MIT-yellow)](LICENSE)
[![Open Access CC BY 4.0](https://img.shields.io/badge/Article-CC%20BY%204.0-lightgrey)](https://creativecommons.org/licenses/by/4.0/)

Data-analysis code and figures for:

> **Longping Fu, Kejun Yin, Xing Xu, and Ronghu Wu.** *Systematic Quantification of Protein O-GlcNAcylation Reveals Common and Cell-Type-Specific Responses to N-Glycosylation Inhibition in Human Cells.* **Analytical Chemistry** 2026, 98, 15689–15699. DOI: [10.1021/acs.analchem.6c00972](https://doi.org/10.1021/acs.analchem.6c00972)

<p align="center">
  <img src="assets/table_of_contents.png" alt="Table of Contents graphic: metabolic labeling and multiplexed proteomics of O-GlcNAcylation in HEK293T, HepG2, and Jurkat cells under N-glycosylation inhibition, resolving common regulation, cell-type-specific responses, and site-specific structural analysis" width="720">
</p>

## Overview

Protein *O*-GlcNAcylation and *N*-glycosylation are two of the most abundant and functionally important glycosylation types in human cells, and both help regulate protein activity, stability, and signaling. *O*-GlcNAcylation — the reversible attachment of a single *N*-acetylglucosamine to serine, threonine, or tyrosine residues of nuclear and cytoplasmic proteins — is a well-established stress sensor, yet how it responds when *N*-glycosylation is perturbed had not been characterized comprehensively or at site-level resolution. Because the two modifications share a nucleotide-sugar precursor (UDP-GlcNAc) and *N*-glycosylation defects trigger ER stress, a crosstalk between them is expected but poorly mapped.

Here, we combined metabolic labeling, bioorthogonal chemistry, and multiplexed (TMT) proteomics to comprehensively and site-specifically quantify protein *O*-GlcNAcylation in three human cell types — **HEK293T, HepG2, and Jurkat** — while inhibiting *N*-glycosylation with **tunicamycin**. Across the cells, more than 1,000 *O*-GlcNAcylated proteins (1,109 total) were identified and quantified. The results reveal both a **common** program shared across cell types and distinct **cell-type-specific** responses to *N*-glycosylation inhibition, and further show that the **local structural context of an *O*-GlcNAc site** is a key determinant of how much it changes. This repository contains the code used to process the mass spectrometry results and generate every figure and supporting table in the paper.

## Highlights

- **Dual-modification-resolved workflow.** Metabolic labeling with Ac₄GalNAz, click-chemistry enrichment, and a galactose-oxidase (GAO) step chemically separate *O*-GlcNAc (+528.2859 Da) from *O*-GalNAc / Tn antigen (+555.2968 Da), a 27 Da mass shift that lets the two structurally similar modifications be distinguished and quantified by MS.
- **1,109 *O*-GlcNAcylated proteins** identified and quantified across HEK293T, HepG2, and Jurkat cells, with 402 commonly detected in all three.
- **Common response to *N*-glycosylation inhibition.** Commonly up-regulated *O*-GlcNAcylated proteins are enriched for translation regulation, glucose response, and RNA splicing; commonly down-regulated ones are enriched for stress-granule assembly and stress-response regulation.
- **Cell-type-specific response.** In Jurkat cells, up-regulated *O*-GlcNAcylated proteins map to leukocyte proliferation, adhesion, and T-cell activation; in HEK293T cells they map to chaperone-mediated folding, cell-cycle, ribonucleotide metabolism, and ribosome biogenesis.
- **Site structure matters.** *O*-GlcNAcylation sites in structured protein regions undergo significantly larger abundance changes than sites in intrinsically disordered regions (IDRs), which stayed largely unaffected — suggesting structured-region sites are more responsive to cellular perturbation.

## Method at a glance

```
Metabolic labeling & treatment   HEK293T / HepG2 / Jurkat  +  tunicamycin (N-glyco inhibition)  +  Ac4GalNAz
        │
Click chemistry & enrichment     CuAAC → Photocleavable-Biotin-Alkyne → NeutrAvidin capture → photocleavage (365 nm)
        │
TMT labeling & GAO oxidation     TMT multiplexing; galactose oxidase distinguishes O-GlcNAc (+528.2859) vs O-GalNAc (+555.2968)
        │
LC-MS/MS                         High-pH HPLC (12 fractions) → Orbitrap Eclipse Tribrid, HCD-pd-EThcD
        │
Database search & analysis       FragPipe + O-Pair Search → differential analysis (|log2(Tuni/Ctrl)| > 0.5, adj. P < 0.05)
```

## Repository structure

```
OGlycoTM/
├── data_analysis/                 R + Python analysis pipeline
│   ├── data_import.R              Import search results
│   ├── data_normalization.R       TMT normalization
│   ├── data_filtering.R           Site / PSM filtering
│   ├── data_quantification.R      O-GlcNAc quantification
│   ├── data_quantification_OGalNAc.R
│   ├── differential_analysis.R    Tuni/Ctrl differential testing
│   ├── glycoprotein_classification.R
│   ├── Figure1.R … Figure6.R      Main-text figure generation
│   ├── FigureS1.R, FigureS2.R     Supporting figures
│   ├── Figure6E_*.py, Figure6F_*.py   PyMOL structure panels (site structural context)
│   ├── generate_supporting_table_S1.R … S11   Supporting tables
│   ├── spectrum_annotator.py,     Glycopeptide spectrum annotation
│   │   fragment_calculator.py
│   └── annotate_*.py, extract_*.py    Spectrum extraction / annotation utilities
├── Manuscript/                    Manuscript, supporting information, revisions
├── assets/                        README figures (TOC graphic)
├── LICENSE                        MIT (code)
└── README.md
```

## Data availability

The mass spectrometry proteomics data have been deposited to the ProteomeXchange Consortium via the **PRIDE** partner repository under the dataset identifier **[PXD073249](https://www.ebi.ac.uk/pride/archive/projects/PXD073249)**.

Per-cell-type identifications and abundance changes for *O*-GlcNAcylated proteins, *O*-GlcNAc sites, *O*-GalNAcylated proteins, and the whole proteome are provided as supporting-information tables with the published article.

## Citation

```bibtex
@article{Fu2026OGlcNAc,
  title   = {Systematic Quantification of Protein O-GlcNAcylation Reveals Common
             and Cell-Type-Specific Responses to N-Glycosylation Inhibition in Human Cells},
  author  = {Fu, Longping and Yin, Kejun and Xu, Xing and Wu, Ronghu},
  journal = {Analytical Chemistry},
  year    = {2026},
  volume  = {98},
  pages   = {15689--15699},
  doi     = {10.1021/acs.analchem.6c00972}
}
```

## Funding

This work was supported by the National Institute of General Medical Sciences of the National Institutes of Health (**R35GM156318**).

## License

The code in this repository is released under the [MIT License](LICENSE). The associated article is open access under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

## Contact

**Ronghu Wu** — School of Chemistry and Biochemistry and the Petit Institute for Bioengineering and Bioscience, Georgia Institute of Technology, Atlanta, GA 30332, USA · [ronghu.wu@chemistry.gatech.edu](mailto:ronghu.wu@chemistry.gatech.edu)
