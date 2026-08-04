from pathlib import Path

import matplotlib.pyplot as plt


OUT_DIR = Path(__file__).resolve().parent

BLACK = "#111111"
BLUE = "#2F6FDB"
ORANGE = "#D97706"
GRAY = "#666666"

ATOM_SIZE = 25
TERMINUS_SIZE = 25
SIDECHAIN_SIZE = 23
ION_SIZE = 27
HYDROGEN_SIZE = 19


def draw_bond(ax, x1, y1, x2, y2, color=BLACK, lw=2.0):
    ax.plot([x1, x2], [y1, y2], color=color, lw=lw, solid_capstyle="round")


def draw_double_bond_vertical(ax, x, y1, y2, color=BLACK, lw=1.7, offset=0.055):
    ax.plot([x - offset, x - offset], [y1, y2], color=color, lw=lw, solid_capstyle="round")
    ax.plot([x + offset, x + offset], [y1, y2], color=color, lw=lw, solid_capstyle="round")


def add_text(ax, x, y, text, size=18, color=BLACK, ha="center", va="center", weight="normal"):
    ax.text(
        x,
        y,
        text,
        ha=ha,
        va=va,
        fontsize=size,
        color=color,
        fontweight=weight,
        family="DejaVu Sans",
    )


def add_arrow(ax, x1, y1, x2, y2, color):
    ax.annotate(
        "",
        xy=(x2, y2),
        xytext=(x1, y1),
        arrowprops=dict(
            arrowstyle="-|>",
            color=color,
            lw=1.6,
            shrinkA=0,
            shrinkB=0,
            mutation_scale=12,
        ),
    )


def dashed_cleavage(ax, x, color, label_top, label_bottom, show_direction=True):
    ax.plot([x, x], [-1.86, 1.86], color=color, lw=1.7, ls=(0, (4, 4)), alpha=0.95)
    add_text(ax, x, 2.20, label_top, size=ION_SIZE, color=color)
    add_text(ax, x, -2.20, label_bottom, size=ION_SIZE, color=color)
    if show_direction:
        # Top ions keep the C terminus; bottom ions keep the N terminus.
        add_arrow(ax, x + 0.17, 2.43, x + 0.52, 2.43, color)
        add_arrow(ax, x - 0.17, -2.43, x - 0.52, -2.43, color)


def draw_backbone(ax):
    """Draw a generic tetrapeptide backbone as a clean presentation schematic."""
    y0 = 0.0
    ca = [1.05, 4.05, 7.05, 10.05]
    carbonyl = [1.95, 4.95, 7.95]
    nitrogens = [2.95, 5.95, 8.95]
    cterm = 11.18

    # N terminus and first alpha carbon.
    add_text(ax, 0.0, y0, r"$\mathrm{H_2N}$", size=TERMINUS_SIZE, ha="right")
    draw_bond(ax, 0.16, y0, ca[0] - 0.28, y0)

    # Residue alpha carbons and generic side chains.
    for i, x in enumerate(ca, start=1):
        add_text(ax, x, y0, "C", size=ATOM_SIZE)
        draw_bond(ax, x, 0.24, x, 0.74)
        draw_bond(ax, x, -0.24, x, -0.70, color=GRAY, lw=1.5)
        add_text(ax, x, 1.06, rf"$\mathrm{{R{i}}}$", size=SIDECHAIN_SIZE)
        add_text(ax, x, -0.98, "H", size=HYDROGEN_SIZE, color=GRAY)

    # Internal carbonyls.
    for x in carbonyl:
        draw_double_bond_vertical(ax, x, 0.26, 0.82)
        add_text(ax, x, 1.03, "O", size=ATOM_SIZE)
        add_text(ax, x, y0, "C", size=ATOM_SIZE)

    # Internal amide nitrogens.
    for x in nitrogens:
        add_text(ax, x, y0, "N", size=ATOM_SIZE)
        draw_bond(ax, x, -0.24, x, -0.70, color=GRAY, lw=1.5)
        add_text(ax, x, -0.98, "H", size=HYDROGEN_SIZE, color=GRAY)

    # Main-chain horizontal bonds.
    atoms = [
        ca[0],
        carbonyl[0],
        nitrogens[0],
        ca[1],
        carbonyl[1],
        nitrogens[1],
        ca[2],
        carbonyl[2],
        nitrogens[2],
        ca[3],
    ]
    for x1, x2 in zip(atoms[:-1], atoms[1:]):
        draw_bond(ax, x1 + 0.28, y0, x2 - 0.28, y0)
    draw_bond(ax, ca[3] + 0.28, y0, cterm - 0.18, y0)

    # C terminus, simplified to keep the cleavage nomenclature visually dominant.
    add_text(ax, cterm, y0, r"$\mathrm{COOH}$", size=TERMINUS_SIZE, ha="left")

    return ca, carbonyl, nitrogens


def make_figure(path, include_cz=False):
    fig, ax = plt.subplots(figsize=(13.4, 5.2), dpi=220)
    fig.patch.set_alpha(0)
    ax.set_facecolor((1, 1, 1, 0))
    ca, carbonyl, nitrogens = draw_backbone(ax)

    # b/y cleavages: peptide amide C-N bond.
    by_positions = [(carbonyl[i] + nitrogens[i]) / 2 for i in range(3)]
    for i, x in enumerate(by_positions, start=1):
        dashed_cleavage(ax, x, BLUE, rf"$y_{{{4 - i}}}$", rf"$b_{{{i}}}$")

    if include_cz:
        # c/z cleavages: N-C_alpha bond after the same residue boundary.
        cz_positions = [(nitrogens[i] + ca[i + 1]) / 2 for i in range(3)]
        for i, x in enumerate(cz_positions, start=1):
            dashed_cleavage(ax, x, ORANGE, rf"$z_{{{4 - i}}}$", rf"$c_{{{i}}}$")

    ax.set_xlim(-0.45, 12.65)
    ax.set_ylim(-2.55, 2.55)
    ax.axis("off")
    fig.savefig(path, transparent=True, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def main():
    make_figure(OUT_DIR / "glycopeptide_fragmentation_by_ions.png", include_cz=False)
    make_figure(OUT_DIR / "glycopeptide_fragmentation_by_cz_ions.png", include_cz=True)


if __name__ == "__main__":
    main()
