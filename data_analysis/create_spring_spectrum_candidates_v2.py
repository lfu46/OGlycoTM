#!/usr/bin/env python3
"""Create checked Spring Symposium spectrum candidate PNGs.

This script fixes HCD/EThcD pairing by verifying the paired HCD precursor
against the EThcD precursor in calibrated mzML metadata. It reads the first
candidate export and writes a new checked folder, leaving the first export
untouched for provenance.
"""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, "/Users/longpingfu/Downloads/GlycoSpectrumAnnotator")
sys.path.insert(0, "/Users/longpingfu/Downloads/mzml-utils/src")

from mzml_utils import MzMLReader  # noqa: E402
from spectrum_annotator_ddzby import (  # noqa: E402
    SpectrumAnnotator,
    parse_modifications_from_string,
)


BASE = Path("/Volumes/cos-lab-rwu60/Longping")
SOURCE = BASE / "OGlycoTM_Final_Version" / "data_source"
SPRING = BASE / "Spring_Symposium_2026"
IN_DIR = SPRING / "Spectrum_Candidates"
OUT_DIR = SPRING / "Spectrum_Candidates_v2"
METRICS_IN = IN_DIR / "candidate_metrics.csv"
FILTER_CACHE = SOURCE / "filter_string_cache.csv"
PPM_TOL = 10.0


def ppm_error(obs: float, ref: float) -> float:
    return abs(obs - ref) / ref * 1e6 if ref else float("inf")


def calibrated_mzml_path(cell_type: str, raw_file: str) -> Path:
    mzml_dir = BASE / "OGlycoTM_Final_Version" / f"OGlycoTM_{cell_type}"
    candidates = [
        mzml_dir / f"{raw_file}_calibrated.mzML",
        mzml_dir / f"{raw_file}_mz_calibrated.mzML",
        mzml_dir / f"{raw_file}.mzML",
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(f"No calibrated mzML found for {raw_file} in {mzml_dir}")


def raw_file_for_ethcd(cache: pd.DataFrame, cell_type: str, ethcd_scan: int) -> str:
    hit = cache[
        (cache["cell_type"] == cell_type)
        & (cache["scan"].astype(int) == int(ethcd_scan))
        & (cache["activation"].astype(str).str.upper() == "ETHCD")
    ]
    if hit.empty:
        raise ValueError(f"No EThcD cache row for {cell_type} scan {ethcd_scan}")
    return str(hit.iloc[0]["raw_file"])


def find_precursor_matched_hcd(reader: MzMLReader, ethcd_scan: int, old_hcd_scan: int) -> tuple[int, str]:
    ethcd = reader.get_spectrum(int(ethcd_scan))
    old = reader.get_spectrum(int(old_hcd_scan))
    if (
        old.activation_type == "HCD"
        and old.precursor_charge == ethcd.precursor_charge
        and ppm_error(float(old.precursor_mz), float(ethcd.precursor_mz)) <= PPM_TOL
    ):
        return int(old_hcd_scan), "kept_existing"

    best_scan = None
    best_key = None
    for scan in range(int(ethcd_scan) - 12, int(ethcd_scan) + 1):
        if scan == int(ethcd_scan):
            continue
        try:
            spec = reader.get_spectrum(scan)
        except Exception:
            continue
        if spec.activation_type != "HCD":
            continue
        if spec.precursor_charge != ethcd.precursor_charge:
            continue
        err = ppm_error(float(spec.precursor_mz), float(ethcd.precursor_mz))
        if err > PPM_TOL:
            continue
        key = (abs(int(ethcd_scan) - scan), err)
        if best_key is None or key < best_key:
            best_scan = scan
            best_key = key
    if best_scan is None:
        raise ValueError(f"No precursor-matched HCD found near EThcD scan {ethcd_scan}")
    return int(best_scan), "corrected_precursor_mismatch"


def glycan_labels_from_mods(mods: list[dict], glycan_text: str) -> dict[int, str]:
    label = "HexNAc-TMT" if "TMT" in str(glycan_text) else "HexNAc"
    labels = {}
    for mod in mods:
        pos = int(mod.get("position", -1))
        mass = float(mod.get("mass", 0.0))
        if pos > 0 and mass > 250:
            labels[pos] = label
    return labels


def ion_label(ion) -> str:
    star = "*" if getattr(ion, "has_modification", False) else ""
    z = "" if ion.charge == 1 else f"^{ion.charge}+"
    loss = f"-{ion.neutral_loss}" if getattr(ion, "neutral_loss", "") else ""
    return f"{ion.ion_type}{ion.ion_number}{star}{loss}{z}"


def bond_coverage(matched, peptide_len: int, ion_types: set[str]) -> tuple[int, str, list[str]]:
    bonds: set[int] = set()
    labels: set[str] = set()
    for ion in matched:
        if ion.ion_type not in ion_types:
            continue
        if getattr(ion, "neutral_loss", ""):
            continue
        if ion.ion_number <= 0:
            continue
        if ion.ion_type in {"b", "c"}:
            bond = int(ion.ion_number)
        else:
            bond = peptide_len - int(ion.ion_number)
        if 1 <= bond < peptide_len:
            bonds.add(bond)
            labels.add(ion_label(ion))
    return len(bonds), f"{len(bonds)}/{peptide_len - 1}", sorted(labels)


def annotate_hcd(row: pd.Series, reader: MzMLReader, hcd_scan: int, mzml_name: str, out_path: Path) -> dict:
    spec = reader.get_spectrum(int(hcd_scan))
    mods = parse_modifications_from_string(str(row["assigned_modifications"]))
    ann = SpectrumAnnotator(
        peptide=str(row["peptide"]),
        modifications=mods,
        precursor_charge=int(spec.precursor_charge),
        precursor_mz=float(spec.precursor_mz),
        exp_mz=spec.mz,
        exp_intensity=spec.intensity,
        tolerance_ppm=20.0,
        site_index=str(row["site_index"]),
        gene=str(row["gene"]),
        activation_type="HCD",
        glycan_labels=glycan_labels_from_mods(mods, str(row.get("glycan", ""))),
        do_deisotope=False,
        scan_num=int(hcd_scan),
        sn_threshold=0.0,
        confidence_level=str(row.get("confidence", "Level1")),
        source_file=mzml_name,
        extended_y_series=False,
    )
    fig = ann.plot(output_path=None)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    found, bonds_text, labels = bond_coverage(ann.matched_ions, int(row["peptide_len"]), {"b", "y"})
    return {
        "hcd_scan": int(hcd_scan),
        "hcd_by_bonds": bonds_text,
        "hcd_by_cov": found / max(int(row["peptide_len"]) - 1, 1),
        "hcd_intensity_annotated": ann.annotation_stats.get("intensity_annotated", ""),
        "hcd_peaks_annotated": ann.annotation_stats.get("peaks_annotated", ""),
        "hcd_fmr_peaks": ann.false_match_rate.fmr_peaks,
        "hcd_n_peaks": len(spec.mz),
        "hcd_by_ions": ",".join(labels),
    }


def copy_single(pattern: str, src_dir: Path, dst_dir: Path) -> None:
    matches = sorted(src_dir.glob(pattern))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected one match for {pattern} in {src_dir}, found {len(matches)}")
    shutil.copy2(matches[0], dst_dir / matches[0].name)


def main() -> None:
    metrics = pd.read_csv(METRICS_IN)
    cache = pd.read_csv(FILTER_CACHE)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    checked_rows = []
    notes = [
        "Spring Symposium 2026 spectrum candidates v2",
        "Generated from Spectrum_Candidates with HCD/EThcD precursor pairing verified by filter_string_cache.csv plus calibrated mzML metadata.",
        "Annotation: spectrum_annotator_ddzby, do_deisotope=False, tolerance=20 ppm, PNG dpi=300.",
        "",
    ]

    readers: dict[Path, MzMLReader] = {}
    try:
        for _, row in metrics.iterrows():
            old_dir = Path(row["candidate_dir"])
            folder_name = old_dir.name
            new_dir = OUT_DIR / folder_name
            new_dir.mkdir(parents=True, exist_ok=True)

            raw_file = raw_file_for_ethcd(cache, str(row["cell_type"]), int(row["ethcd_scan"]))
            mzml_path = calibrated_mzml_path(str(row["cell_type"]), raw_file)
            if mzml_path not in readers:
                readers[mzml_path] = MzMLReader(str(mzml_path))
            reader = readers[mzml_path]

            checked_hcd, status = find_precursor_matched_hcd(
                reader, int(row["ethcd_scan"]), int(row["hcd_scan"])
            )

            copy_single("*_MS1_isolation.png", old_dir, new_dir)
            copy_single(f"*_EThcD_s{int(row['ethcd_scan'])}.png", old_dir, new_dir)

            new_row = row.copy()
            new_row["candidate_dir"] = str(new_dir)
            if checked_hcd == int(row["hcd_scan"]):
                copy_single(f"*_HCD_s{int(row['hcd_scan'])}.png", old_dir, new_dir)
            else:
                hcd_out = new_dir / f"{folder_name}_HCD_s{checked_hcd}.png"
                updates = annotate_hcd(row, reader, checked_hcd, mzml_path.name, hcd_out)
                for key, value in updates.items():
                    new_row[key] = value

            checked_rows.append(new_row)
            notes.append(
                f"{row['gene']} {row['site_index']} ({row['cell_type']}, EThcD s{int(row['ethcd_scan'])}): "
                f"HCD s{int(row['hcd_scan'])} -> s{checked_hcd} ({status}); dir {new_dir}"
            )
    finally:
        for reader in readers.values():
            reader.close()

    checked = pd.DataFrame(checked_rows)
    checked.to_csv(OUT_DIR / "candidate_metrics_checked.csv", index=False)
    (OUT_DIR / "candidate_notes_checked.txt").write_text("\n".join(notes) + "\n")
    print(f"Wrote checked candidates to {OUT_DIR}")
    print(checked[["gene", "site_index", "cell_type", "ethcd_scan", "hcd_scan", "hcd_by_bonds", "ethcd_cz_bonds", "oxonium_max_rel"]].to_string(index=False))


if __name__ == "__main__":
    main()
