# PyMOL 3D Structure Visualization (Figure 5E/5F, formerly 6E/6F)

> Moved out of `CLAUDE.md` on 2026-09-02 so it is loaded on demand rather than into every session.

Site panels are **table-driven**. The 16 `OGalNAc_*_pymol.py` and 4 `Figure6E_*_pymol.py` scripts
were two templates repeated with three values changed; they were replaced 2026-08-04 by
`pymol_site_panels.py` + `pymol_site_panels.csv` (2,143 lines → 309), verified to render
pixel-identically.

```bash
python3 pymol_site_panels.py --list          # show the panel table
python3 pymol_site_panels.py --only PTPRC    # render one (matches gene/accession/site)
python3 pymol_site_panels.py                 # render all 20
```

To add a panel, add a CSV row — do not write a new script. Columns:
`style,gene,accession,sites,color,out_subdir,prefix,notes`, where `style` is `surface`
(whole residue as spheres, the O-GalNAc panels) or `sidechain` (only the modified side-chain
atoms, the Figure 6E panels), and `sites` is a semicolon list like `T139;S146`.

Shared formatting (unchanged from the originals):
- pLDDT coloring (blue=high confidence, orange=low)
- Transparent surface (70% transparency)
- Site highlighted as colored spheres (orange=upregulated, cyan=stable)

**NEVER hand-write an AlphaFold filename or version.** The driver resolves structures via
`mzml_utils.structure.fetch_structure()` with a cached-file fast path. Interpolating
`AF-{acc}-F1-model_v6.pdb` is how `Figure6F_EWSR1_S274_pymol.py` ended up pinned to `model_v4`
while every other panel used v6 — that script now resolves the newest cached model instead.

Still hand-written, kept for genuine per-protein domain/IDR colouring that does not parameterize:
`Figure6F_EWSR1_S274_pymol.py`, `Figure6F_HOXA13_pymol.py`, `Figure6F_HYOU1_pymol.py`,
`Figure6F_pymol.py`, and `Figure6F_ray_modes.py` (render-settings exploration).

```bash
# Run a hand-written PyMOL script (requires PyMOL installed via Homebrew)
/opt/homebrew/bin/pymol -c -q Figure6F_HYOU1_pymol.py
```

PyMOL ray trace modes:
- Mode 0: Default (clean, no shadows)
- Mode 1: Shadows (realistic, publication quality)
- Mode 2: Black outline only (line drawing)
- Mode 3: Quantized/posterized colors
