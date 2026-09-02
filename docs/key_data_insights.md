# Key Data Insights

> Moved out of `CLAUDE.md` on 2026-09-02 so it is loaded on demand rather than into every session.

- **O-GlcNAc sites**: ~90% are in IDR (intrinsically disordered regions), ~10% in structured regions
- **IDR classification**: Uses StructureMap's pPSE method - smoothed `nAA_24_180_pae` (pPSE_24_smooth10) ≤ 34.27 = IDR. This is calculated in `structuremap_analysis.py` following the StructureMap tutorial. Note: pLDDT is included as a regression predictor but is NOT used for IDR classification.
- **Site features file**: `site_features/OGlcNAc_site_features.csv` contains per-site logFC, pLDDT, is_IDR, secondary structure, pPSE_24, pPSE_12, pPSE_24_smooth10
- **Proteins with sites in both regions**: Very rare (only 3-4 proteins across all cell types)
- **OGlcNAc Atlas reference**: `reference/OGlcNAcAtlas_unambiguous_sites_20251208.csv` for checking reported vs novel sites
