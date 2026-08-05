#!/usr/bin/env python3
"""
export_web_data.py — build the compact JSON bundle that powers the OGlycoTM
interactive data browser (docs/ static site on GitHub Pages).

Reads the processed differential-expression, normalization, and site-feature
CSVs from the network drive and emits small, rounded JSON files into
docs/data/. The data is derived/publication-level (the same content as the
Supporting Tables and the PRIDE deposition), so it is safe to commit.

Run:  python3 data_analysis/export_web_data.py
"""
import csv
import json
import math
import os
import re

# ---------------------------------------------------------------- config
# The network share (SMB) drops reads intermittently; stage CSVs locally and
# point OGLYCO_DATA at the copy. Falls back to the share if unset.
BASE = os.environ.get(
    "OGLYCO_DATA",
    "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source",
)
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "docs", "data")
CELLS = ["HEK293T", "HepG2", "Jurkat"]

LFC_THRESH = 0.5        # |log2(Tuni/Ctrl)| threshold used in the paper
ADJP_THRESH = 0.05      # Benjamini-Hochberg adjusted p threshold

CHANNELS = ["Tuni_1", "Tuni_2", "Tuni_3", "Ctrl_4", "Ctrl_5", "Ctrl_6"]

# Row counts of the published Supporting Tables, in CELLS order. The browser is
# a second public view of the same results, so it has to agree with them exactly;
# main() fails rather than publish a bundle that does not.
PUBLISHED_COUNTS = {
    "oglcnac_protein": ("S2", [773, 673, 671]),
    "oglcnac_site": ("S6", [313, 257, 226]),
    "ogalnac_protein": ("S8", [275, 173, 245]),
    "ogalnac_site": ("S10", [171, 80, 154]),
    "wp_protein": ("S4", [9201, 8104, 8672]),
}

os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------- helpers
def read_csv(path):
    if not os.path.exists(path):
        print(f"  !! missing: {path}")
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def fnum(x):
    """Parse to float or None (NaN/blank/inf -> None)."""
    if x is None:
        return None
    s = str(x).strip()
    if s == "" or s.lower() in ("na", "nan", "inf", "-inf"):
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    if math.isnan(v) or math.isinf(v):
        return None
    return v


def sig(x, n=3):
    """Round to n significant figures, keeping None."""
    v = fnum(x)
    if v is None:
        return None
    if v == 0:
        return 0.0
    return float(f"{v:.{n}g}")


def rnd(x, d=3):
    v = fnum(x)
    return None if v is None else round(v, d)


def sig_flag(lfc, adjp):
    """1 = up, -1 = down, 0 = not significant."""
    lfc, adjp = fnum(lfc), fnum(adjp)
    if lfc is None or adjp is None:
        return 0
    if adjp < ADJP_THRESH and abs(lfc) > LFC_THRESH:
        return 1 if lfc > 0 else -1
    return 0


def pick_intensity_cols(fieldnames):
    """Return the 6 intensity column names in CHANNELS order, preferring the
    fully-normalized (_sl_tmm) columns, then _sl, then raw."""
    for suffix in ("_sl_tmm", "_sl", ""):
        cols = [f"Intensity.{ch}{suffix}" for ch in CHANNELS]
        if all(c in fieldnames for c in cols):
            return cols, (suffix or "raw")
    return None, None


def intensities_lookup(path, key_col):
    """Map key -> [6 rounded normalized intensities] from a norm CSV."""
    rows = read_csv(path)
    if not rows:
        return {}, None
    cols, kind = pick_intensity_cols(rows[0].keys())
    if cols is None:
        return {}, None
    out = {}
    for r in rows:
        vals = [sig(r.get(c), 4) for c in cols]
        out[r.get(key_col)] = vals
    return out, kind


SITE_RE = re.compile(r"^(?P<prot>[A-Z0-9\-]+)_(?P<res>[A-Z])(?P<pos>\d+)")


def parse_site_index(site_index):
    """P35658_T1363 -> ('P35658', 'T', 1363)."""
    m = SITE_RE.match(str(site_index))
    if not m:
        return None, None, None
    return m.group("prot"), m.group("res"), int(m.group("pos"))


# ------------------------------------------------------- site confidence filter
# The published site tables keep Level1, plus Level1b whose O-Pair localization
# probability is >= 0.75 (fix_all_tables.R:351 for S5/S6, update_supporting_
# tables.R:287 for S9/S10). The site DE CSVs carry neither column, so the
# accepted set is rebuilt here from the site PSM tables -- without it the browser
# publishes low-confidence sites the paper excluded.
SITE_PROB_MIN = 0.75
PROB_RE = re.compile(r"[0-9.]+(?=\]$)")


def parse_site_probability(prob_str):
    """Trailing localization probability of an O-Pair Site.Probabilities string."""
    m = PROB_RE.search(str(prob_str or ""))
    if not m:
        return None
    try:
        return float(m.group(0))
    except ValueError:
        return None


def high_conf_site_index(glyco):
    """cell -> set of site_index values that pass the published filter."""
    allow = {}
    for cell in CELLS:
        keep = set()
        for r in read_csv(os.path.join(BASE, "site", f"{glyco}_site_{cell}.csv")):
            level = (r.get("Confidence.Level") or "").strip()
            if level == "Level1":
                keep.add(r.get("site_index"))
            elif level == "Level1b":
                p = parse_site_probability(r.get("Site.Probabilities"))
                if p is not None and p >= SITE_PROB_MIN:
                    keep.add(r.get("site_index"))
        allow[cell] = keep
    return allow


def write_json(name, obj):
    path = os.path.join(OUT, name)
    with open(path, "w") as f:
        json.dump(obj, f, separators=(",", ":"), allow_nan=False)
    kb = os.path.getsize(path) / 1024
    print(f"  wrote {name:<32} {kb:8.1f} KB")
    return kb


# ---------------------------------------------------------------- gene index
# gene -> layer -> cell -> summary, for the global search / per-gene card
gene_index = {}


def gi(gene):
    g = (gene or "").strip()
    if not g:
        g = "(unnamed)"
    return gene_index.setdefault(g, {})


def set_meta(gene, pid, name, desc):
    """Record a display name/desc/id for a gene once (first, highest-priority
    layer wins — protein layers are built before sites/WP)."""
    m = gi(gene)
    if "m" not in m:
        m["m"] = [pid or "", name or "", desc or ""]


# ---------------------------------------------------------------- protein layers
def build_protein_layer(tag, de_prefix, norm_prefix, gi_key):
    """O-GlcNAc / O-GalNAc / WP protein DE -> long-format records + gene index."""
    records = []
    for cell in CELLS:
        de = read_csv(os.path.join(BASE, "differential_analysis", f"{de_prefix}_{cell}.csv"))
        norm, kind = intensities_lookup(
            os.path.join(BASE, "normalization", f"{norm_prefix}_{cell}.csv"),
            "UniProt_Accession" if de_prefix.startswith("WP") else "Protein.ID",
        )
        for r in de:
            # WP rows carry a clean UniProt_Accession alongside a "sp|ACC|NAME"
            # Protein.ID; prefer the clean accession so links resolve.
            pid = r.get("UniProt_Accession") or r.get("Protein.ID")
            gene = r.get("Gene") or r.get("Gene.Symbol") or ""
            lfc = rnd(r.get("logFC"), 3)
            adjp = sig(r.get("adj.P.Val"), 3)
            pval = sig(r.get("P.Value"), 3)
            flag = sig_flag(r.get("logFC"), r.get("adj.P.Val"))
            key = r.get("UniProt_Accession") if de_prefix.startswith("WP") else pid
            rec = {
                "id": pid,
                "gene": gene,
                "name": r.get("Entry.Name") or "",
                "desc": (r.get("Protein.Description") or r.get("Annotation") or "").split(" OS=")[0].strip(),
                "cell": cell,
                "lfc": lfc,
                "adjP": adjp,
                "p": pval,
                "sig": flag,
                "int": norm.get(key),
            }
            records.append(rec)
            g = gi(gene).setdefault(gi_key, {})
            g[cell] = {"id": pid, "lfc": lfc, "adjP": adjp, "sig": flag}
            set_meta(gene, pid, rec["name"], rec["desc"])
    return records


# ---------------------------------------------------------------- O-GlcNAc sites
def build_oglcnac_sites():
    allow = high_conf_site_index("OGlcNAc")
    feats = read_csv(os.path.join(BASE, "site_features", "OGlcNAc_site_features.csv"))
    # extra descriptive fields (Entry.Name, Protein.Description) from site DE
    de_extra = {}
    for cell in CELLS:
        for r in read_csv(os.path.join(BASE, "differential_analysis", f"OGlcNAc_site_DE_{cell}.csv")):
            de_extra[r.get("site_index")] = {
                "name": r.get("Entry.Name") or "",
                "desc": r.get("Protein.Description") or "",
            }
    norm = {}
    for cell in CELLS:
        lk, _ = intensities_lookup(
            os.path.join(BASE, "normalization", f"OGlcNAc_site_norm_{cell}.csv"), "site_index"
        )
        norm.update(lk)

    records = []
    for r in feats:
        si = r.get("site_index")
        cell = r.get("cell")
        if si not in allow.get(cell, ()):
            continue
        prot, res, pos = parse_site_index(si)
        if res is None:
            res = "S" if r.get("is_serine") in ("1", "1.0") else "T"
            pos = int(fnum(r.get("site_number")) or 0)
        lfc = rnd(r.get("logFC"), 3)
        adjp = sig(r.get("adj.P.Val"), 3)
        flag = sig_flag(r.get("logFC"), r.get("adj.P.Val"))
        extra = de_extra.get(si, {})
        rec = {
            "site": si,
            "id": r.get("Protein.ID"),
            "gene": r.get("Gene") or "",
            "name": extra.get("name", ""),
            "desc": extra.get("desc", ""),
            "cell": cell,
            "res": res,
            "pos": pos,
            "peptide": r.get("Peptide") or "",
            "lfc": lfc,
            "adjP": adjp,
            "sig": flag,
            "plddt": rnd(r.get("pLDDT"), 1),
            "idr": 1 if str(r.get("is_IDR")).startswith("1") else 0,
            "ss": r.get("secondary_structure_simple") or "",
            "ppse": rnd(r.get("pPSE_24_smooth10"), 2),
            "hyd": rnd(r.get("hydrophobicity_7mer"), 2),
            "int": norm.get(si),
        }
        records.append(rec)
        g = gi(rec["gene"]).setdefault("og_s", [])
        g.append({"site": si, "cell": cell, "lfc": lfc, "adjP": adjp,
                  "idr": rec["idr"], "res": res, "pos": pos})
        set_meta(rec["gene"], rec["id"], rec["name"], rec["desc"])
    return records


# ---------------------------------------------------------------- O-GalNAc sites
def build_ogalnac_sites():
    allow = high_conf_site_index("OGalNAc")
    norm = {}
    for cell in CELLS:
        lk, _ = intensities_lookup(
            os.path.join(BASE, "normalization", f"OGalNAc_site_norm_{cell}.csv"), "site_index"
        )
        norm.update(lk)
    records = []
    for cell in CELLS:
        for r in read_csv(os.path.join(BASE, "differential_analysis", f"OGalNAc_site_DE_{cell}.csv")):
            si = r.get("site_index")
            if si not in allow.get(cell, ()):
                continue
            prot, res, pos = parse_site_index(si)
            if res is None:
                res = (r.get("modified_residue") or "?")[:1]
                pos = int(fnum(r.get("site_number")) or 0)
            lfc = rnd(r.get("logFC"), 3)
            adjp = sig(r.get("adj.P.Val"), 3)
            flag = sig_flag(r.get("logFC"), r.get("adj.P.Val"))
            gene = r.get("Gene") or ""
            # intensities: prefer norm lookup, else inline _sl columns
            ints = norm.get(si)
            if ints is None:
                cols, _ = pick_intensity_cols(r.keys())
                if cols:
                    ints = [sig(r.get(c), 4) for c in cols]
            rec = {
                "site": si, "id": r.get("Protein.ID"), "gene": gene,
                "desc": r.get("Protein.Description") or "", "cell": cell,
                "res": res, "pos": pos, "lfc": lfc, "adjP": adjp, "sig": flag,
                "int": ints,
            }
            records.append(rec)
            g = gi(gene).setdefault("ga_s", [])
            g.append({"site": si, "cell": cell, "lfc": lfc, "adjP": adjp, "res": res, "pos": pos})
            set_meta(gene, r.get("Protein.ID"), "", rec["desc"])
    return records


# ---------------------------------------------------------------- main
def counts_by_cell(records):
    out = {c: {"n": 0, "up": 0, "down": 0} for c in CELLS}
    for r in records:
        c = r["cell"]
        out[c]["n"] += 1
        if r.get("sig") == 1:
            out[c]["up"] += 1
        elif r.get("sig") == -1:
            out[c]["down"] += 1
    return out


def check_against_published(counts):
    """Abort if any layer disagrees with its Supporting Table row count."""
    bad = []
    for layer, (table, expected) in PUBLISHED_COUNTS.items():
        got = [counts[layer][c]["n"] for c in CELLS]
        mark = "ok" if got == expected else "MISMATCH"
        print(f"  {layer:<18} {str(got):<22} Table {table:<3} {str(expected):<22} {mark}")
        if got != expected:
            bad.append(f"{layer}: browser {got} vs Table {table} {expected}")
    if bad:
        raise SystemExit(
            "\nRefusing to write a bundle that disagrees with the published tables:\n  "
            + "\n  ".join(bad)
        )


def main():
    print("Building web data bundle -> docs/data/")

    print("O-GlcNAc proteins ...")
    og_p = build_protein_layer("O-GlcNAc protein", "OGlcNAc_protein_DE",
                               "OGlcNAc_protein_norm", "og_p")
    print("O-GlcNAc sites ...")
    og_s = build_oglcnac_sites()
    print("O-GalNAc proteins ...")
    ga_p = build_protein_layer("O-GalNAc protein", "OGalNAc_protein_DE",
                               "OGalNAc_protein_norm", "ga_p")
    print("O-GalNAc sites ...")
    ga_s = build_ogalnac_sites()
    print("Whole proteome ...")
    wp = build_protein_layer("WP protein", "WP_protein_DE",
                             "WP_protein_norm", "wp")

    meta = {
        "cells": CELLS,
        "lfc_thresh": LFC_THRESH,
        "adjp_thresh": ADJP_THRESH,
        "counts": {
            "oglcnac_protein": counts_by_cell(og_p),
            "oglcnac_site": counts_by_cell(og_s),
            "ogalnac_protein": counts_by_cell(ga_p),
            "ogalnac_site": counts_by_cell(ga_s),
            "wp_protein": counts_by_cell(wp),
        },
        "n_genes": len(gene_index),
        # union across cell types -- the overview stat strip reads this rather
        # than carrying its own copy of the number
        "n_oglcnac_sites_union": len({r["site"] for r in og_s}),
    }

    print("Checking against the published Supporting Tables ...")
    check_against_published(meta["counts"])

    sizes = {}
    for name, payload in (
        ("oglcnac_protein.json", og_p),
        ("oglcnac_site.json", og_s),
        ("ogalnac_protein.json", ga_p),
        ("ogalnac_site.json", ga_s),
        ("wp_protein.json", wp),
        ("gene_index.json", gene_index),
    ):
        sizes[name] = write_json(name, payload)
    write_json("meta.json", meta)

    print(f"\nTotal bundle: {sum(sizes.values())/1024:.2f} MB across {len(sizes)+1} files")
    print(f"Unique genes indexed: {len(gene_index)}")
    print("Counts (n / up / down) by cell:")
    for layer, cc in meta["counts"].items():
        s = "  ".join(f"{c}={v['n']}/{v['up']}/{v['down']}" for c, v in cc.items())
        print(f"  {layer:<18} {s}")


if __name__ == "__main__":
    main()
