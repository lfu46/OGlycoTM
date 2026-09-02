# Storage layout: file map, artifact regeneration and drive sync

> Moved out of `CLAUDE.md` on 2026-09-02 so it is loaded on demand rather than into every session.

Read **`00_FILE_MAP.md` at the network project root** before touching the data. It records the
layout, the search provenance (FragPipe 24.0 / MSFragger 4.4 / pwiz 3.0.25139), and the exact
commands that regenerate the two deleted classes of artifact:

- 43 GB of 600-DPI TIFF spectra — rebuilt from the kept PDFs with `magick -density 600 … -quality 100`
  (verified pixel-identical).
- 44.7 GB of **uncalibrated** `.mzML` — rebuilt from the kept `.raw` with the recorded `msconvert`
  filter chain. Analysis uses the *calibrated* mzML, which was kept.

Drift between the two drives is declared in `oglycotm_sync_policy.tsv` on the Expansion drive and
checked with the shared checker:

```bash
python3 /Users/longpingfu/Downloads/OGlyco_DBA/data_analysis/tools/dba_sync_check.py \
        --policy /Volumes/Expansion/Longping/OGlycoTM/oglycotm_sync_policy.tsv
#   ... --checksum   md5-verify same-size files (catches SMB truncation)
#   ... --apply      reconcile (never deletes; a size mismatch needs a human)
```
