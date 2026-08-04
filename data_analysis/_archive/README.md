# Archived scripts

Superseded or exhausted scripts from the Anal. Chem. revision round, moved here 2026-08-04.
Nothing is deleted — these are the record of how the published analysis was arrived at, and
`git log --follow` still works through the move.

**None of these are part of the reproducible pipeline.** For that, see the scripts still in
`data_analysis/`.

## Superseded by a newer script in the same folder

| Archived | Superseded by | Why |
|---|---|---|
| `generate_supporting_table_S1.R` … `S7.R` | `fix_all_tables.R` | Rebuilds S1–S7 from source with correct numeric types into `version_Apr5_2026/`; the old ones wrote the stale `new_version/` under pre-revision numbering |
| `generate_supporting_table_S8.R` | `regenerate_S11.R` | S-GlcNAc/Cys content moved from slot S8 to S11; the old script would now overwrite S8 (O-GalNAc protein DE) with the wrong content |
| `generate_supporting_table_S9_OGalNAc.R` | `update_supporting_tables.R` (Part 3) | Writes the identical file |
| `generate_supporting_table_S10_OGalNAc_DE.R` | `update_supporting_tables.R` (Part 2) | O-GalNAc protein DE is S8 in the revised numbering, not S10 |
| `Figure6.R` | `Figure6_filtered.R` | Same panels without the Level1b site-probability ≥ 0.75 filter; its 6E/6F examples were replaced |
| `OGalNAc_site_ranking_plot.R` | `OGalNAc_site_norm_and_plot.R` | Ranking plot on raw rather than SL+TMM-normalized intensities |
| `OGalNAc_plots.R` | `OGalNAc_site_norm_and_plot.R` + `OGalNAc_GSEA_selected_plot.R` | Both of its two parts have dedicated newer replacements |
| `OGalNAc_secretory_down_barplots.R` | `Figure6D_secretory_barplot.R` | Exploratory 5-category version replaced by the 9-protein / 3-category panel that shipped |
| `find_ideal_examples.R` | `find_ideal_sites.R` | Protein-average search replaced by a site-level search across all three cell types |

`Figure5.R` is a special case: the subcellular-localization panel it builds was **removed from
the manuscript** in the March 2026 revision (old Figure 6 became the new Figure 5). It is kept
because its `wilcox_across_locations()` / `wilcox_across_cells()` / `prepare_stat_annotations()`
helpers — auto-stacked significance brackets with BH correction and `min_n` guards — are the best
statistics scaffolding in the repo and are flagged for extraction into `wu-lab-r`.

## Superseded by an installed library

| Archived | Use instead |
|---|---|
| `add_activation_type.py`, `extract_activation_type.py` | `mzml_utils.classify_activation` + `parse_spectrum_id` + `find_mzml_file` (also handles the MS1 case, which neither local copy did) |
| `extract_spectra.py`, `extract_ethcd_spectra.py` | `mzml_utils.open_spectra()` / `MzMLReader.get_spectrum()` |
| `annotate_low_prob_spectra.py` | `annotate_low_prob_hcd_spectra.py` — written 49 minutes later against the installed `spectrum_annotator_ddzby` |
| `reannotate_tyr_spectra.py`, `generate_tyr_spectra_v3.py` | `regenerate_tyr_spectra_v5.py` |

The 20 hand-written `*_pymol.py` site scripts were not archived but **replaced** —
see `pymol_site_panels.py` + `pymol_site_panels.csv`, verified to render pixel-identically.

## Exhausted one-offs

`check_site_filter.R`, `check_outliers.R`, `check_reference_formats.R`, `check_site_examples.R`,
`verify_all_numbers.R`, `analyze_low_prob_psm.R`, `find_all_both.R`, `find_alternative.R`,
`find_ideal_sites.R`, `compare_filtering_approaches.R`, `OGalNAc_GSEA_exclusive.R`,
`final_check.R` — each answered a single question during the revision; the answers are baked into
the shipped scripts and tables.

`final_check.R` is **not runnable today** regardless: it reads `/tmp/ms_final.txt`,
`/tmp/p2p_final.txt` and `/tmp/si_final.txt`, which no longer exist.
