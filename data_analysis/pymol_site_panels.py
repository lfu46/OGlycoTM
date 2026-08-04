#!/usr/bin/env python3
"""Render O-GalNAc / O-GlcNAc site structure panels from a table.

Replaces 20 near-identical hand-written `*_pymol.py` scripts (~1,750 lines) that differed only
in accession, gene name, residue and output path. The panels they produce are unchanged: the
PyMOL command sequences below are transcribed from the two templates those scripts shared, so
re-rendering reproduces the published figures rather than restyling them.

Two things are deliberately NOT carried over from the originals:

  * The AlphaFold path is no longer hand-written. Every original interpolated
    `AF-{acc}-F1-model_v6.pdb` -- and Figure6F_EWSR1_S274 still said `model_v4`, which no longer
    exists upstream. Structures now come from `mzml_utils.structure.fetch_structure()`, which
    resolves the correct current file through 3D-Beacons and caches it.
  * Output directory resolution (env var -> network -> local fallback) applied to all panels,
    not just the two that happened to have it.

Usage:
    python3 pymol_site_panels.py --list            # show the table
    python3 pymol_site_panels.py --dry-run         # resolve structures, emit no images
    python3 pymol_site_panels.py                   # render everything
    python3 pymol_site_panels.py --only PTPRC      # substring match on gene/accession/sites

Env:
    OGLYCOTM_FIGURES    Figures root (default: the network share)
    OGLYCOTM_STRUCTURES structure cache (default: <figures-root>/../data_source/alphafold_structures)
"""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
TABLE = HERE / "pymol_site_panels.csv"

NET_FIGURES = Path("/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures")
NET_STRUCTURES = Path(
    "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/alphafold_structures"
)
FALLBACK_FIGURES = HERE / "Figure_panels_preview"

# Side-chain atoms carrying the hydroxyl that gets glycosylated, by residue letter.
SIDECHAIN_ATOMS = {
    "S": ["CB", "OG"],
    "T": ["CB", "OG1", "CG2"],
    "Y": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
}

PLDDT_SETUP = """
cmd.set_color("plddt_very_high", [0/255, 83/255, 214/255])
cmd.set_color("plddt_high", [101/255, 203/255, 243/255])
cmd.set_color("plddt_low", [255/255, 219/255, 19/255])
cmd.set_color("plddt_very_low", [255/255, 125/255, 69/255])
stored.plddt_colors = {{}}
cmd.iterate("{obj} and name CA", "stored.plddt_colors[resi] = b")
for resi, plddt in stored.plddt_colors.items():
    if plddt > 90: cmd.color("plddt_very_high", f"{obj} and resi {{resi}}")
    elif plddt > 70: cmd.color("plddt_high", f"{obj} and resi {{resi}}")
    elif plddt > 50: cmd.color("plddt_low", f"{obj} and resi {{resi}}")
    else: cmd.color("plddt_very_low", f"{obj} and resi {{resi}}")
"""


def resolve_dir(env_var: str, default: Path, fallback: Path) -> Path:
    """env var -> default (if writable) -> local fallback. Lifted from the two scripts that had it."""
    env = os.environ.get(env_var)
    if env:
        p = Path(env)
        p.mkdir(parents=True, exist_ok=True)
        return p
    try:
        default.mkdir(parents=True, exist_ok=True)
        if os.access(default, os.W_OK):
            return default
    except OSError:
        pass
    fallback.mkdir(parents=True, exist_ok=True)
    return fallback


def load_table(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        lines = [ln for ln in fh if not ln.lstrip().startswith("#")]
    return [r for r in csv.DictReader(lines) if r.get("gene")]


def get_structure(accession: str, cache: Path) -> Path | None:
    """Cached local file first (keeps this runnable offline), else the structure resolver.

    Never builds an AlphaFold URL by hand -- fetch_structure() goes through 3D-Beacons and
    returns the correct current version, which interpolating `model_v6` cannot guarantee.
    """
    if cache.is_dir():
        hits = sorted(cache.glob(f"AF-{accession}-F1-model_v*.pdb"), reverse=True)
        if hits:
            return hits[0]
    try:
        from mzml_utils.structure import fetch_structure
    except ImportError:
        print(f"    mzml_utils.structure unavailable and no cached file for {accession}")
        return None
    try:
        res = fetch_structure(accession, cache_dir=cache)
    except Exception as exc:                                    # network/resolver failure
        print(f"    fetch_structure({accession}) failed: {exc}")
        return None
    path = getattr(res, "path", None)
    if path is None:
        print(f"    no structure available for {accession} (status={getattr(res, 'status', '?')})")
        return None
    return Path(path)


def build_script(row: dict, pdb: Path, out_dir: Path) -> tuple[str, str]:
    obj = row["gene"]
    sites = [s.strip() for s in row["sites"].split(";") if s.strip()]
    stem = f"{row['prefix']}_{obj}_{'_'.join(sites)}"
    color = row["color"]
    body = [
        "from pymol import cmd, stored",
        "import os, subprocess",
        f'cmd.load(r"{pdb}", "{obj}")',
        'cmd.bg_color("white")',
        'cmd.set("ray_opaque_background", 1)',
        'cmd.set("antialias", 2)',
        'cmd.set("ray_trace_mode", 1)',
        'cmd.set("ray_shadows", 0)',
        'cmd.set("depth_cue", 0)',
        'cmd.set("fog", 0)',
    ]

    if row["style"] == "surface":
        body += ['cmd.set("spec_reflect", 0.3)', 'cmd.set("spec_power", 200)']
    body += ['cmd.hide("everything")']
    body.append(PLDDT_SETUP.format(obj=obj))
    body += [f'cmd.show("cartoon", "{obj}")']
    if row["style"] == "surface":
        body += ['cmd.set("cartoon_fancy_helices", 1)', 'cmd.set("cartoon_smooth_loops", 1)']
    body += [
        f'cmd.show("surface", "{obj}")',
        f'cmd.set("transparency", 0.7, "{obj}")',
        'cmd.set("surface_quality", 1)',
    ]

    for site in sites:
        aa, pos = site[0], site[1:]
        sel = f"site_{site}"
        body.append(f'cmd.select("{sel}", "{obj} and resi {pos}")')
        if row["style"] == "surface":
            body += [
                f'cmd.show("spheres", "{sel}")',
                f'cmd.color("{color}", "{sel}")',
                f'cmd.set("sphere_scale", 2.0, "{sel}")',
                f'cmd.set("sphere_transparency", 0, "{sel}")',
            ]
        else:
            atoms = " or ".join(f"name {a}" for a in SIDECHAIN_ATOMS.get(aa, ["CB"]))
            body += [
                f'cmd.show("spheres", "{sel} and ({atoms})")',
                f'cmd.color("{color}", "{sel}")',
                f'cmd.set("sphere_scale", 0.8, "{sel}")',
            ]
    body.append("cmd.deselect()")

    if row["style"] == "surface":
        body += [
            f'cmd.orient("{obj}")', f'cmd.zoom("{obj}", buffer=5)',
            'cmd.turn("y", 45)', 'cmd.turn("x", -20)',
            'cmd.set("ambient", 0.3)', 'cmd.set("direct", 0.7)',
            'cmd.set("specular", 0.4)', 'cmd.set("shininess", 40)',
            'cmd.set("ray_shadows", 1)',
            'cmd.set("ray_shadow_decay_factor", 0.3)',
            'cmd.set("ray_shadow_decay_range", 2.0)',
            'cmd.set("light_count", 2)', 'cmd.set("light", "[-0.4, -0.4, -1.0]")',
        ]
    else:
        body += [
            f'cmd.orient("{obj}")', 'cmd.turn("y", 30)', 'cmd.turn("x", -10)',
            f'cmd.zoom("{obj}", buffer=5)',
            'cmd.set("ambient", 0.35)', 'cmd.set("direct", 0.6)',
            'cmd.set("specular", 0.3)', 'cmd.set("ray_shadows", 1)',
        ]

    png = out_dir / f"{stem}_structure.png"
    pse = out_dir / f"{stem}_session.pse"
    pdf = out_dir / f"{stem}_structure.pdf"
    body += [
        'cmd.set("antialias", 2)',
        f'cmd.save(r"{pse}")',
        "cmd.ray(1200, 1200)",
        f'cmd.png(r"{png}", dpi=600)',
    ]
    if row["style"] == "surface":
        # sips is macOS-only; the originals used it too, so failure is non-fatal.
        body.append(
            f'subprocess.run(["sips", "-s", "format", "pdf", r"{png}", "--out", r"{pdf}"],'
            " capture_output=True)"
        )
    body.append(f'print("wrote {stem}")')
    return "\n".join(body) + "\n", stem


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", type=Path, default=TABLE)
    ap.add_argument("--only", help="substring match on gene, accession or sites")
    ap.add_argument("--list", action="store_true", help="print the table and exit")
    ap.add_argument("--dry-run", action="store_true", help="resolve structures, render nothing")
    args = ap.parse_args()

    rows = load_table(args.table)
    if args.only:
        q = args.only.lower()
        rows = [r for r in rows
                if q in r["gene"].lower() or q in r["accession"].lower() or q in r["sites"].lower()]
        if not rows:
            print(f"no row matches {args.only!r}", file=sys.stderr)
            return 2

    if args.list:
        for r in rows:
            print(f"  {r['style']:10} {r['gene']:10} {r['accession']:8} {r['sites']:18} "
                  f"{r['color']:10} {r['out_subdir']}")
        print(f"\n{len(rows)} panel(s)")
        return 0

    figures = resolve_dir("OGLYCOTM_FIGURES", NET_FIGURES, FALLBACK_FIGURES)
    structures = Path(os.environ.get("OGLYCOTM_STRUCTURES", NET_STRUCTURES))
    print(f"figures    {figures}\nstructures {structures}\n")

    try:
        from mzml_utils.structure.render import pymol_binary
        pymol = pymol_binary()
    except Exception:
        pymol = "/opt/homebrew/bin/pymol"

    ok = failed = 0
    for row in rows:
        label = f"{row['gene']} {row['sites']}"
        print(f"[{label}]")
        pdb = get_structure(row["accession"], structures)
        if pdb is None or not Path(pdb).is_file():
            failed += 1
            continue
        out_dir = figures / row["out_subdir"]
        out_dir.mkdir(parents=True, exist_ok=True)
        script, stem = build_script(row, Path(pdb), out_dir)
        if args.dry_run:
            print(f"    would render {stem} from {Path(pdb).name}")
            ok += 1
            continue
        with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as fh:
            fh.write(script)
            tmp = fh.name
        try:
            proc = subprocess.run([pymol, "-c", "-q", tmp], capture_output=True, text=True)
            if proc.returncode != 0:
                print(f"    pymol exit={proc.returncode}: {proc.stderr.strip()[:300]}")
                failed += 1
            else:
                print(f"    {stem}_structure.png")
                ok += 1
        finally:
            os.unlink(tmp)

    print(f"\n{ok} ok, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
