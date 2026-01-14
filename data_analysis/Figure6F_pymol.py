"""
Figure 6F - TAB2 (Q9NYJ8) 3D Structure with O-GlcNAc Sites
PyMOL Python script for visualizing AlphaFold structure colored by pLDDT
with O-GlcNAc sites highlighted as spheres

Run this script in PyMOL:
    pymol Figure6F_pymol.py

Or from within PyMOL:
    run Figure6F_pymol.py
"""

from pymol import cmd, stored
import os

# ============================================
# Configuration
# ============================================

# File paths
STRUCTURE_PATH = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/alphafold_structures/AF-Q9NYJ8-F1-model_v6.pdb"
OUTPUT_DIR = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/Figure6"

# O-GlcNAc site data (HEK293T)
SITES = {
    29: {"logFC": 0.824, "region": "Structured", "pLDDT": 88.2},
    359: {"logFC": -0.272, "region": "IDR", "pLDDT": 31.8},
    460: {"logFC": -1.54, "region": "IDR", "pLDDT": 53.3}
}

# Colors
COLOR_UP = "tv_orange"      # For positive logFC (up-regulated)
COLOR_DOWN = "cyan"         # For negative logFC (down-regulated)
COLOR_LABEL = "black"

# ============================================
# Main Script
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
    cmd.load(STRUCTURE_PATH, "TAB2")
    cmd.hide("everything")
    cmd.show("cartoon", "TAB2")

def color_by_plddt():
    """Color structure by pLDDT confidence score (B-factor)"""
    # AlphaFold pLDDT scheme:
    # Very high (>90): dark blue #0053D6
    # High (70-90): light blue #65CBF3
    # Low (50-70): yellow #FFDB13
    # Very low (<50): orange #FF7D45

    # Define custom colors (RGB values 0-1 scale)
    cmd.set_color("plddt_very_high", [0/255, 83/255, 214/255])      # >90
    cmd.set_color("plddt_high", [101/255, 203/255, 243/255])        # 70-90
    cmd.set_color("plddt_low", [255/255, 219/255, 19/255])          # 50-70
    cmd.set_color("plddt_very_low", [255/255, 125/255, 69/255])     # <50

    # Use alter to color by B-factor ranges
    # PyMOL's spectrum doesn't give the exact AlphaFold colors, so we iterate
    stored.plddt_colors = {}

    # Get B-factors and assign colors
    cmd.iterate("TAB2 and name CA", "stored.plddt_colors[resi] = b")

    # Color each residue based on pLDDT
    for resi, plddt in stored.plddt_colors.items():
        if plddt > 90:
            cmd.color("plddt_very_high", f"TAB2 and resi {resi}")
        elif plddt > 70:
            cmd.color("plddt_high", f"TAB2 and resi {resi}")
        elif plddt > 50:
            cmd.color("plddt_low", f"TAB2 and resi {resi}")
        else:
            cmd.color("plddt_very_low", f"TAB2 and resi {resi}")

def highlight_oglcnac_sites():
    """Highlight O-GlcNAc sites as spheres"""
    for resi, data in SITES.items():
        site_name = f"site_S{resi}"

        # Select the serine residue
        cmd.select(site_name, f"TAB2 and resi {resi}")

        # Show as spheres (side chain atoms)
        cmd.show("spheres", f"{site_name} and (name CB or name OG)")

        # Color based on fold change direction
        if data["logFC"] > 0:
            cmd.color(COLOR_UP, f"{site_name} and (name CB or name OG)")
        else:
            cmd.color(COLOR_DOWN, f"{site_name} and (name CB or name OG)")

        # Set sphere size
        cmd.set("sphere_scale", 1.2, site_name)

    # Create combined selection
    cmd.select("all_sites", " or ".join([f"site_S{r}" for r in SITES.keys()]))

def add_site_labels():
    """Add labels showing site position and fold change"""
    cmd.set("label_font_id", 7)  # Sans-serif bold
    cmd.set("label_size", 20)
    cmd.set("label_color", COLOR_LABEL)
    cmd.set("label_outline_color", "white")

    for resi, data in SITES.items():
        fc_str = f"{data['logFC']:.2f}"
        if data['logFC'] > 0:
            fc_str = f"+{fc_str}"

        # Label at CA atom position
        label_text = f"S{resi} (FC={fc_str})"
        cmd.label(f"TAB2 and resi {resi} and name CA", f'"{label_text}"')

def setup_view():
    """Set up the camera view"""
    cmd.orient("TAB2")
    cmd.turn("y", 45)
    cmd.turn("x", -15)
    cmd.zoom("TAB2", buffer=8)

def render_and_save():
    """Ray trace and save the image in multiple formats"""
    # ============================================
    # Mode 1: Shadows - Realistic with soft shadows
    # ============================================

    cmd.set("ray_trace_mode", 1)

    # Lighting settings
    cmd.set("ambient", 0.35)
    cmd.set("direct", 0.6)
    cmd.set("specular", 0.3)
    cmd.set("ray_shadows", 1)
    cmd.set("ray_shadow_decay_factor", 0.1)

    # Anti-aliasing
    cmd.set("antialias", 2)
    cmd.set("ray_trace_frames", 1)
    cmd.set("ray_opaque_background", 1)

    # Depth cueing off
    cmd.set("depth_cue", 0)
    cmd.set("fog", 0)

    # Save session for manual adjustment
    session_path = os.path.join(OUTPUT_DIR, "Figure6F_TAB2_session.pse")
    cmd.save(session_path)
    print(f"Session saved: {session_path}")

    # Output size: 2 inch x 2 inch at 600 DPI = 1200 x 1200 pixels (high resolution)
    width_px = 1200
    height_px = 1200
    output_dpi = 600

    # === Version 1: With labels (PNG, ray-traced) ===
    cmd.ray(width_px, height_px)
    png_with_labels = os.path.join(OUTPUT_DIR, "Figure6F_TAB2_with_labels.png")
    cmd.png(png_with_labels, dpi=output_dpi)
    print(f"PNG with labels saved: {png_with_labels}")

    # === Version 2: Without labels (PNG, ray-traced) ===
    # Hide labels temporarily
    cmd.hide("labels")
    cmd.ray(width_px, height_px)
    png_no_labels = os.path.join(OUTPUT_DIR, "Figure6F_TAB2_no_labels.png")
    cmd.png(png_no_labels, dpi=output_dpi)
    print(f"PNG without labels saved: {png_no_labels}")

    # === Version 3: EPS (vector format, editable text) ===
    # Note: EPS doesn't support ray tracing, but text is editable
    cmd.show("labels")  # Re-show labels for EPS
    eps_path = os.path.join(OUTPUT_DIR, "Figure6F_TAB2_structure.eps")
    cmd.set("ray_default_renderer", 0)  # Disable ray for vector output
    cmd.png(eps_path.replace(".eps", "_preview.png"), width_px, height_px, dpi=output_dpi)  # Preview
    # Use multisave for EPS
    try:
        cmd.multisave(eps_path, "TAB2", format="eps")
        print(f"EPS (vector) saved: {eps_path}")
    except:
        print("EPS export not available in this PyMOL version")

    # === Version 4: Convert PNG to PDF ===
    pdf_path = os.path.join(OUTPUT_DIR, "Figure6F_TAB2_structure.pdf")
    try:
        import subprocess
        # Use ImageMagick convert or sips (macOS built-in)
        # First try sips (native macOS)
        result = subprocess.run(
            ["sips", "-s", "format", "pdf", png_no_labels, "--out", pdf_path],
            capture_output=True, text=True
        )
        if result.returncode == 0:
            print(f"PDF (from PNG) saved: {pdf_path}")
        else:
            # Try convert (ImageMagick)
            result = subprocess.run(
                ["convert", png_no_labels, pdf_path],
                capture_output=True, text=True
            )
            if result.returncode == 0:
                print(f"PDF saved: {pdf_path}")
            else:
                print("PDF conversion tools not available")
    except Exception as e:
        print(f"PDF conversion not available: {e}")

    print("\n** For best results with editable text: **")
    print("1. Open the 'no_labels' PNG in Adobe Illustrator or Inkscape")
    print("2. Add text labels manually as vector text")
    print("3. Save as PDF")

def create_legend_info():
    """Print legend information for reference"""
    print("\n" + "="*50)
    print("Figure 6F - TAB2 O-GlcNAc Sites")
    print("="*50)
    print("\nModel confidence (pLDDT):")
    print("  - Very high (>90): Dark blue")
    print("  - High (70-90): Light blue")
    print("  - Low (50-70): Yellow")
    print("  - Very low (<50): Orange")
    print("\nO-GlcNAc sites (spheres):")
    print(f"  - Orange: Up-regulated (logFC > 0)")
    print(f"  - Cyan: Down-regulated (logFC < 0)")
    print("\nSite details:")
    for resi, data in SITES.items():
        direction = "UP" if data["logFC"] > 0 else "DOWN"
        print(f"  S{resi}: logFC={data['logFC']:.2f} ({direction}), "
              f"Region={data['region']}, pLDDT={data['pLDDT']}")
    print("="*50 + "\n")

# ============================================
# Execute
# ============================================

if __name__ == "pymol" or __name__ == "__main__":
    print("Generating Figure 6F...")

    setup_pymol()
    load_structure()
    color_by_plddt()
    highlight_oglcnac_sites()
    add_site_labels()
    setup_view()
    create_legend_info()

    # Auto-render the image
    render_and_save()

    print("\nStructure loaded! Adjust view manually if needed, then run:")
    print("  render_and_save()")
    print("\nOr render manually:")
    print("  ray 2400, 2400")
    print("  png /Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/Figure6F_TAB2_structure.png, dpi=300")
