"""
Figure 6E - DDX50 Y492 (Q9BQ39) 3D Structure - HepG2
O-GlcNAc site in Helicase C-terminal domain (380-524)

Site: Y492
  - logFC: +0.95 (significantly upregulated)
  - pLDDT: 94.94 (very high confidence - structured)
  - Domain: Helicase C-terminal (380-524)
  - Cell type: HepG2

Run: /opt/homebrew/bin/pymol -c -q Figure6E_DDX50_Y492_pymol.py
"""

from pymol import cmd, stored
import os

# ============================================
# Configuration
# ============================================

PROTEIN_ID = "Q9BQ39"
PROTEIN_NAME = "DDX50"
STRUCTURE_PATH = f"/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/alphafold_structures/AF-{PROTEIN_ID}-F1-model_v6.pdb"
OUTPUT_DIR = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/Figure6/Figure6E_candidates"

# O-GlcNAc site data
SITE_POSITION = 492
SITE_AA = "Y"
SITE_DATA = {
    "logFC": 0.95,
    "pLDDT": 94.94,
    "domain": "Helicase C-terminal",
    "domain_range": (380, 524),
    "cell_type": "HepG2"
}

# Colors
SITE_COLOR = "tv_orange"  # Orange for upregulated site

# ============================================
# Functions
# ============================================

def setup_pymol():
    """Initialize PyMOL settings"""
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_shadows", 0)
    cmd.set("depth_cue", 0)
    cmd.set("fog", 0)

def load_structure():
    """Load and prepare the AlphaFold structure"""
    cmd.load(STRUCTURE_PATH, PROTEIN_NAME)
    cmd.hide("everything")
    cmd.show("cartoon", PROTEIN_NAME)

def color_by_plddt():
    """Color structure by pLDDT confidence score"""
    cmd.set_color("plddt_very_high", [0/255, 83/255, 214/255])
    cmd.set_color("plddt_high", [101/255, 203/255, 243/255])
    cmd.set_color("plddt_low", [255/255, 219/255, 19/255])
    cmd.set_color("plddt_very_low", [255/255, 125/255, 69/255])

    stored.plddt_colors = {}
    cmd.iterate(f"{PROTEIN_NAME} and name CA", "stored.plddt_colors[resi] = b")

    for resi, plddt in stored.plddt_colors.items():
        if plddt > 90:
            cmd.color("plddt_very_high", f"{PROTEIN_NAME} and resi {resi}")
        elif plddt > 70:
            cmd.color("plddt_high", f"{PROTEIN_NAME} and resi {resi}")
        elif plddt > 50:
            cmd.color("plddt_low", f"{PROTEIN_NAME} and resi {resi}")
        else:
            cmd.color("plddt_very_low", f"{PROTEIN_NAME} and resi {resi}")

def add_surface():
    """Add transparent surface representation with pLDDT coloring"""
    cmd.show("surface", PROTEIN_NAME)
    cmd.set("transparency", 0.7, PROTEIN_NAME)
    cmd.set("surface_quality", 1)

def highlight_domain():
    """Highlight the domain containing the O-GlcNAc site"""
    domain_start, domain_end = SITE_DATA["domain_range"]
    cmd.select("domain_region", f"{PROTEIN_NAME} and resi {domain_start}-{domain_end}")
    cmd.deselect()

def highlight_oglcnac_site():
    """Highlight O-GlcNAc site as sphere"""
    site_name = f"site_{SITE_AA}{SITE_POSITION}"

    # Select the residue
    cmd.select(site_name, f"{PROTEIN_NAME} and resi {SITE_POSITION}")

    # Show side chain as spheres (Tyrosine - phenol ring)
    cmd.show("spheres", f"{site_name} and (name CB or name CG or name CD1 or name CD2 or name CE1 or name CE2 or name CZ or name OH)")

    # Color the site
    cmd.color(SITE_COLOR, site_name)
    cmd.set("sphere_scale", 0.8, site_name)

    cmd.deselect()

def setup_view():
    """Set up the camera view - focus on the site"""
    cmd.orient(PROTEIN_NAME)
    # Rotate to show the Helicase C-terminal domain clearly
    cmd.turn("y", 45)
    cmd.turn("x", -15)
    cmd.zoom(PROTEIN_NAME, buffer=5)

def render_and_save():
    """Ray trace and save the image"""
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    cmd.set("ray_trace_mode", 1)
    cmd.set("ambient", 0.35)
    cmd.set("direct", 0.6)
    cmd.set("specular", 0.3)
    cmd.set("ray_shadows", 1)
    cmd.set("antialias", 2)

    width_px = 1200
    height_px = 1200
    output_dpi = 600

    session_path = os.path.join(OUTPUT_DIR, f"Figure6E_{PROTEIN_NAME}_{SITE_AA}{SITE_POSITION}_session.pse")
    cmd.save(session_path)
    print(f"Session saved: {session_path}")

    cmd.ray(width_px, height_px)
    png_path = os.path.join(OUTPUT_DIR, f"Figure6E_{PROTEIN_NAME}_{SITE_AA}{SITE_POSITION}_structure.png")
    cmd.png(png_path, dpi=output_dpi)
    print(f"PNG saved: {png_path}")

def print_info():
    """Print site information"""
    print("\n" + "="*50)
    print(f"Figure 6E - {PROTEIN_NAME} {SITE_AA}{SITE_POSITION} ({SITE_DATA['cell_type']})")
    print("="*50)
    print(f"\nSite: {SITE_AA}{SITE_POSITION}")
    print(f"  logFC: +{SITE_DATA['logFC']:.2f} (upregulated)")
    print(f"  pLDDT: {SITE_DATA['pLDDT']} (structured)")
    print(f"  Domain: {SITE_DATA['domain']} ({SITE_DATA['domain_range'][0]}-{SITE_DATA['domain_range'][1]})")
    print("\nColor scheme:")
    print("  - Structure: pLDDT confidence (blue=high, orange=low)")
    print("  - Surface: 70% transparent, pLDDT coloring")
    print(f"  - Site {SITE_AA}{SITE_POSITION}: Orange sphere (upregulated)")
    print("="*50 + "\n")

# ============================================
# Execute
# ============================================

if __name__ == "pymol" or __name__ == "__main__":
    print(f"Generating Figure 6E - {PROTEIN_NAME} {SITE_AA}{SITE_POSITION}...")

    setup_pymol()
    load_structure()
    color_by_plddt()
    add_surface()
    highlight_domain()
    highlight_oglcnac_site()
    setup_view()
    print_info()
    render_and_save()

    print("\nDone! Structure rendered.")
