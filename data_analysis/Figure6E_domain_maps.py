#!/usr/bin/env python3
"""
Figure 6E: Protein domain maps with O-GalNAc site annotations.
Inspired by GlycoDomainViewer (Joshi et al., Glycobiology 2018).
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
matplotlib.rcParams['font.family'] = 'Arial'

DEFAULT_OUTPUT_DIR = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/Figure6_OGalNAc"
FALLBACK_OUTPUT_DIR = "/Users/longpingfu/Downloads/OGlycoTM/data_analysis/Figure6_OGalNAc_preview"


def resolve_output_dir():
    env_output_dir = os.environ.get("O_GALNAC_OUTPUT_DIR")
    if env_output_dir:
        os.makedirs(env_output_dir, exist_ok=True)
        return env_output_dir

    os.makedirs(DEFAULT_OUTPUT_DIR, exist_ok=True)
    if os.access(DEFAULT_OUTPUT_DIR, os.W_OK):
        return DEFAULT_OUTPUT_DIR

    os.makedirs(FALLBACK_OUTPUT_DIR, exist_ok=True)
    return FALLBACK_OUTPUT_DIR


OUTPUT_DIR = resolve_output_dir()

# Color scheme
# Project cell-type colors: HEK293T=#4DBBD5, HepG2=#F39B7F, Jurkat=#00A087
DOMAIN_COLORS = {
    "Mucin-like": "#A8DAE5",             # HEK293T light
    "Cys-rich": "#4DBBD5",               # HEK293T
    "Fibronectin type-III": "#F39B7F",   # HepG2
    "Tyrosine-protein phosphatase": "#00A087",  # Jurkat
    "Peptidase M14": "#4DBBD5",          # HEK293T
    "MRH": "#00A087",                    # Jurkat
}
TM_COLOR = "#FFD700"
SP_COLOR = "#B0BEC5"
BACKBONE_COLOR = "#E0E0E0"
SITE_FILL_COLOR = "#FFD700"
SITE_EDGE_COLOR = "#222222"
SITE_TEXT_COLOR = "#E64B35"
NOTE_COLOR = "#5F6368"

proteins = [
    {
        "gene": "CPD",
        "uniprot": "O75976",
        "length": 1380,
        "length_label_offset": 45,
        "signal": (1, 31),
        "show_signal_label": True,
        "tm": (1300, 1320),
        "domains": [
            (57, 380, "Peptidase\nM14-1", "Peptidase M14"),
            (502, 792, "Peptidase\nM14-2", "Peptidase M14"),
            (932, 1211, "Peptidase\nM14-3", "Peptidase M14"),
        ],
        "site_region_note": None,
        "note_fontsize": 7.8,
        "note_y_offset": -0.55,
        "topology_labels": [
            (630, "Luminal domain", "#4477AA", "center", -0.41),
            (1415, "Cytosolic tail", "#AA4444", "left", -0.41),
        ],
        "sites": [
            {
                "pos": 44,
                "label": "T44",
                "ratio": 0.09,
                "label_side": "right",
                "marker_x_offset": 0.0,
                "marker_y_offset": 0.30,
                "text_y_offset": 0.30,
                "text_x_offset": 0.02,
                "placement": "top",
            },
        ],
    },
    {
        "gene": "PTPRC (CD45)",
        "uniprot": "P08575",
        "length": 1306,
        "length_label_offset": 45,
        "signal": (1, 25),
        "show_signal_label": False,
        "tm": (578, 598),
        "domains": [
            (26, 227, "Mucin", "Mucin-like"),
            (228, 390, "Cys", "Cys-rich"),
            (391, 483, "FN1", "Fibronectin type-III"),
            (484, 576, "FN2", "Fibronectin type-III"),
            (653, 912, "PTP D1", "Tyrosine-protein phosphatase"),
            (944, 1228, "PTP D2", "Tyrosine-protein phosphatase"),
        ],
        "topology_labels": [
            (255, "Extracellular", "#4477AA", "center", -0.52),
            (950, "Cytoplasmic", "#AA4444", "center", -0.52),
        ],
        "sites": [
            {
                "pos": 139,
                "label": "T139",
                "ratio": 0.59,
                "label_side": "left",
                "marker_x_offset": 0.0,
                "marker_y_offset": 0.30,
                "text_y_offset": 0.30,
                "text_x_offset": -0.04,
                "placement": "top",
            },
            {
                "pos": 146,
                "label": "S146",
                "ratio": 0.27,
                "label_side": "right",
                "marker_x_offset": 0.0,
                "marker_y_offset": -0.30,
                "text_y_offset": -0.30,
                "text_x_offset": 0.04,
                "placement": "bottom",
            },
        ],
    },
]


def add_box(ax, x, y, width, height, facecolor, edgecolor, linewidth, zorder):
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle="round,pad=0.003",
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=linewidth,
            zorder=zorder,
        )
    )


def draw_protein(ax, prot, y_center, max_len=None, x_offset=0.12, usable_width=0.86, bar_height=0.30):
    """Draw one protein domain map. Each protein fills the same width."""
    length = prot["length"]
    scale = usable_width / length  # each protein fills full width
    bar_start = x_offset
    bar_width = usable_width
    backbone_y = y_center - bar_height / 2

    # Backbone.
    add_box(
        ax,
        bar_start,
        backbone_y,
        bar_width,
        bar_height,
        BACKBONE_COLOR,
        "#9A9A9A",
        0.8,
        1,
    )

    # Signal peptide.
    if prot.get("signal"):
        s, e = prot["signal"]
        add_box(
            ax,
            bar_start + s * scale,
            backbone_y,
            (e - s) * scale,
            bar_height,
            SP_COLOR,
            "#8C8C8C",
            0.6,
            2,
        )
        if prot.get("show_signal_label", True):
            ax.text(
                bar_start + ((s + e) / 2) * scale,
                y_center,
                "SP",
                ha="center",
                va="center",
                fontsize=4,
                fontweight="bold",
                color="white",
                zorder=3,
            )

    # Domains.
    for s, e, label, dtype in prot["domains"]:
        domain_width = (e - s) * scale
        add_box(
            ax,
            bar_start + s * scale,
            backbone_y,
            (e - s) * scale,
            bar_height,
            DOMAIN_COLORS.get(dtype, "#BDBDBD"),
            "#575757",
            0.8,
            4,
        )
        ax.text(
            bar_start + ((s + e) / 2) * scale,
            y_center,
            label,
            ha="center",
            va="center",
            fontsize=5,
            fontweight="bold",
            color="black",
            zorder=5,
        )

    # Transmembrane segment.
    if prot.get("tm"):
        s, e = prot["tm"]
        tm_height = bar_height * 1.5
        tm_y = y_center - tm_height / 2
        add_box(
            ax,
            bar_start + s * scale,
            tm_y,
            (e - s) * scale,
            tm_height,
            TM_COLOR,
            "#B8860B",
            0.8,
            6,
        )
        # Use same y_offset as topology labels for consistent alignment
        tm_y_offset = prot["topology_labels"][0][4] if len(prot["topology_labels"][0]) == 5 else prot.get("topology_y_offset", -0.31)
        ax.text(
            bar_start + ((s + e) / 2) * scale,
            y_center + tm_y_offset,
            "TM",
            ha="center",
            va="top",
            fontsize=6,
            color="#555555",
            style="italic",
            zorder=7,
        )

    # Site annotations.
    for site in prot["sites"]:
        x = bar_start + site["pos"] * scale
        marker_x = x + site.get("marker_x_offset", 0.0)
        marker_y = y_center + site["marker_y_offset"]
        text_y = y_center + site["text_y_offset"]
        placement = site.get("placement", "top")
        stem_start_y = y_center + bar_height / 2 if placement == "top" else y_center - bar_height / 2
        stem_end_y = marker_y - 0.025 if placement == "top" else marker_y + 0.025
        text_x = marker_x + site.get(
            "text_x_offset",
            -0.02 if site["label_side"] == "left" else 0.02,
        )

        ax.plot(
            [x, x],
            [stem_start_y, stem_end_y],
            color="#444444",
            linewidth=0.8,
            zorder=8,
        )
        ax.plot(
            x,
            marker_y,
            marker="s",
            markersize=5,
            color=SITE_FILL_COLOR,
            markeredgecolor=SITE_EDGE_COLOR,
            markeredgewidth=0.8,
            zorder=9,
        )

        text = f"{site['label']} ({site['ratio']:.2f})"
        if site["label_side"] == "left":
            ax.text(
                text_x,
                text_y,
                text,
                ha="right",
                va="center",
                fontsize=6,
                fontweight="bold",
                color=SITE_TEXT_COLOR,
                zorder=10,
            )
        else:
            ax.text(
                text_x,
                text_y,
                text,
                ha="left",
                va="center",
                fontsize=6,
                fontweight="bold",
                color=SITE_TEXT_COLOR,
                zorder=10,
            )

    # Protein name and length.
    ax.text(
        x_offset - 0.03,
        y_center,
        prot["gene"],
        ha="right",
        va="center",
        fontsize=6,
        fontweight="bold",
        color="black",
    )
    # aa number labels at bottom of each protein
    ax.text(bar_start, y_center - bar_height / 2 - 0.03, '1',
            ha='center', va='top', fontsize=4.5, color='#666666')
    ax.text(bar_start + usable_width, y_center - bar_height / 2 - 0.03, str(length),
            ha='center', va='top', fontsize=4.5, color='#666666')

    # Topology labels.
    for topology in prot["topology_labels"]:
        if len(topology) == 5:
            pos, label, color, ha, y_offset = topology
        else:
            pos, label, color, ha = topology
            y_offset = prot.get("topology_y_offset", -0.31)
        ax.text(
            bar_start + pos * scale,
            y_center + y_offset,
            label,
            ha=ha,
            va="top",
            fontsize=6,
            color=color,
            style="italic",
        )

    # Site-region context note.
    if prot.get("site_region_note"):
        site_midpoint = prot.get(
            "note_pos",
            sum(site["pos"] for site in prot["sites"]) / len(prot["sites"]),
        )
        ax.text(
            bar_start + site_midpoint * scale,
            y_center + prot.get("note_y_offset", -0.49),
            prot["site_region_note"],
            ha=prot.get("note_ha", "center"),
            va="center",
            fontsize=prot.get("note_fontsize", 8.5),
            color=NOTE_COLOR,
        )


def main():
    fig, ax = plt.subplots(figsize=(3.8, 1.5))
    ax.set_xlim(0, 1.08)
    ax.set_ylim(-0.4, 2.1)
    ax.axis("off")

    y_positions = [1.45, 0.42]
    max_len = max(prot["length"] for prot in proteins)
    for prot, y in zip(proteins, y_positions):
        draw_protein(ax, prot, y, max_len)

    # GalNAc legend — yellow square + label
    legend_y = -0.1
    legend_x_text = 1.08  # right-align with "Cytosolic tail"
    legend_x_marker = legend_x_text - 0.12
    ax.plot(legend_x_marker, legend_y, 's', markersize=5, color='#FFD700',
            markeredgecolor='black', markeredgewidth=0.8, zorder=10)
    ax.text(legend_x_text, legend_y, 'GalNAc', ha='right', va='center',
            fontsize=6, color='#333333', zorder=10)
    pass  # no border

    fig.tight_layout()
    pdf_path = os.path.join(OUTPUT_DIR, "Figure6E_domain_maps.pdf")
    png_path = os.path.join(OUTPUT_DIR, "Figure6E_domain_maps.png")
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved Figure6E to {pdf_path}")
    print(f"Saved Figure6E to {png_path}")


if __name__ == "__main__":
    main()
