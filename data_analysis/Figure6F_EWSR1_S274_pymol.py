"""
Figure 6F - EWSR1 S274 (Q01844) 3D Structure - HEK293T
O-GlcNAc site in IDR region

Site: S274
  - logFC: -0.15 (minimal change - stable)
  - pLDDT: 44.5 (low confidence - IDR)
  - Region: Intrinsically Disordered Region (IDR)
  - Cell type: HEK293T

Run: /opt/homebrew/bin/pymol -c -q Figure6F_EWSR1_S274_pymol.py
"""

from pymol import cmd, stored
import os

# ============================================
# Configuration
# ============================================

PROTEIN_ID = "Q01844"
PROTEIN_NAME = "EWSR1"
STRUCTURE_PATH = f"/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/alphafold_structures/AF-{PROTEIN_ID}-F1-model_v4.pdb"
OUTPUT_DIR = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/Figure6"

# O-GlcNAc site data
SITE_POSITION = 274
SITE_AA = "S"
SITE_DATA = {
    "logFC": -0.15,
    "pLDDT": 44.5,
    "region": "IDR",
    "cell_type": "HEK293T"
}

# Colors
SITE_COLOR = "cyan"  # Cyan for IDR site with minimal change

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
    """Add transparent surface representation"""
    cmd.show("surface", PROTEIN_NAME)
    cmd.set("transparency", 0.7, PROTEIN_NAME)
    cmd.set("surface_quality", 1)

def highlight_oglcnac_site():
    """Highlight O-GlcNAc site as sphere"""
    site_name = f"site_{SITE_AA}{SITE_POSITION}"

    # Select the residue
    cmd.select(site_name, f"{PROTEIN_NAME} and resi {SITE_POSITION}")

    # Show side chain as spheres (Serine)
    cmd.show("spheres", f"{site_name} and (name CB or name OG)")

    # Color the site
    cmd.color(SITE_COLOR, site_name)
    cmd.set("sphere_scale", 1.5, site_name)
    cmd.set("sphere_transparency", 0, site_name)  # Site sphere is opaque

    cmd.deselect()

def setup_view():
    """Set up the camera view - focus on the site"""
    cmd.orient(PROTEIN_NAME)
    cmd.turn("y", 30)
    cmd.turn("x", -10)
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

    session_path = os.path.join(OUTPUT_DIR, f"Figure6F_{PROTEIN_NAME}_{SITE_AA}{SITE_POSITION}_session.pse")
    cmd.save(session_path)
    print(f"Session saved: {session_path}")

    cmd.ray(width_px, height_px)
    png_path = os.path.join(OUTPUT_DIR, f"Figure6F_{PROTEIN_NAME}_{SITE_AA}{SITE_POSITION}_structure.png")
    cmd.png(png_path, dpi=output_dpi)
    print(f"PNG saved: {png_path}")

def print_info():
    """Print site information"""
    print("\n" + "="*50)
    print(f"Figure 6F - {PROTEIN_NAME} {SITE_AA}{SITE_POSITION} ({SITE_DATA['cell_type']})")
    print("="*50)
    print(f"\nSite: {SITE_AA}{SITE_POSITION}")
    print(f"  logFC: {SITE_DATA['logFC']:.2f} (minimal change)")
    print(f"  pLDDT: {SITE_DATA['pLDDT']} (IDR)")
    print(f"  Region: {SITE_DATA['region']}")
    print("\nColor scheme:")
    print("  - Structure: pLDDT confidence (blue=high, orange=low)")
    print("  - Surface: 70% transparent")
    print(f"  - Site {SITE_AA}{SITE_POSITION}: Cyan sphere (IDR site)")
    print("="*50 + "\n")

# ============================================
# Execute
# ============================================

if __name__ == "pymol" or __name__ == "__main__":
    print(f"Generating Figure 6F - {PROTEIN_NAME} {SITE_AA}{SITE_POSITION}...")

    setup_pymol()
    load_structure()
    color_by_plddt()
    add_surface()
    highlight_oglcnac_site()
    setup_view()
    print_info()
    render_and_save()

    print("\nDone! Structure rendered.")
