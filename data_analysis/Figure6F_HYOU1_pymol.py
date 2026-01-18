"""
Figure 6F - HYOU1 (Q9Y4L1) 3D Structure with O-GlcNAc Site T871
PyMOL Python script for visualizing AlphaFold structure with transparent surface
showing secondary structure underneath, colored by pLDDT

Run this script in PyMOL:
    pymol Figure6F_HYOU1_pymol.py

Or from command line:
    /opt/homebrew/bin/pymol -c -q Figure6F_HYOU1_pymol.py
"""

from pymol import cmd, stored
import os

# ============================================
# Configuration
# ============================================

STRUCTURE_PATH = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/alphafold_structures/AF-Q9Y4L1-F1-model_v6.pdb"
OUTPUT_DIR = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/Figure6"

# O-GlcNAc site data (HEK293T) - STRUCTURED REGION EXAMPLE
# T871: Threonine at position 871, in structured region (pLDDT = 94.1)
SITE = {
    "residue": 871,
    "residue_type": "T",  # Threonine
    "logFC": 2.47,
    "pLDDT": 94.1,
    "region": "Structured"
}

# Colors
COLOR_SITE = "tv_orange"  # Orange for the O-GlcNAc site (up-regulated)
SURFACE_TRANSPARENCY = 0.7  # 0 = opaque, 1 = fully transparent

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
    cmd.set("spec_reflect", 0.3)
    cmd.set("spec_power", 200)

def load_structure():
    """Load the AlphaFold structure"""
    cmd.load(STRUCTURE_PATH, "HYOU1")
    cmd.hide("everything")

def color_by_plddt():
    """Color structure by pLDDT confidence score (stored in B-factor)"""
    # AlphaFold pLDDT color scheme
    cmd.set_color("plddt_very_high", [0/255, 83/255, 214/255])      # >90: Dark blue
    cmd.set_color("plddt_high", [101/255, 203/255, 243/255])        # 70-90: Light blue
    cmd.set_color("plddt_low", [255/255, 219/255, 19/255])          # 50-70: Yellow
    cmd.set_color("plddt_very_low", [255/255, 125/255, 69/255])     # <50: Orange

    # Get B-factors (pLDDT) for each residue
    stored.plddt_colors = {}
    cmd.iterate("HYOU1 and name CA", "stored.plddt_colors[resi] = b")

    # Color each residue based on pLDDT
    for resi, plddt in stored.plddt_colors.items():
        if plddt > 90:
            cmd.color("plddt_very_high", f"HYOU1 and resi {resi}")
        elif plddt > 70:
            cmd.color("plddt_high", f"HYOU1 and resi {resi}")
        elif plddt > 50:
            cmd.color("plddt_low", f"HYOU1 and resi {resi}")
        else:
            cmd.color("plddt_very_low", f"HYOU1 and resi {resi}")

def show_cartoon_and_surface():
    """Show cartoon with transparent surface overlay"""
    # Show cartoon (secondary structure)
    cmd.show("cartoon", "HYOU1")
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_smooth_loops", 1)

    # Show transparent surface
    cmd.show("surface", "HYOU1")
    cmd.set("transparency", SURFACE_TRANSPARENCY, "HYOU1")
    cmd.set("surface_quality", 1)

def highlight_oglcnac_site():
    """Highlight O-GlcNAc site T871 as large spheres"""
    resi = SITE["residue"]

    # Select the threonine residue - show entire residue as spheres
    cmd.select("site_T871", f"HYOU1 and resi {resi}")

    # Show the entire residue as spheres (all atoms)
    cmd.show("spheres", f"site_T871")

    # Color the site
    cmd.color(COLOR_SITE, f"site_T871")

    # Set larger sphere size for visibility
    cmd.set("sphere_scale", 2.0, "site_T871")

    # Make site visible through surface (no transparency on spheres)
    cmd.set("sphere_transparency", 0, "site_T871")

def setup_view():
    """Set up the camera view - show overall structure"""
    # Orient to show entire protein
    cmd.orient("HYOU1")

    # Zoom to show the complete structure with some buffer
    cmd.zoom("HYOU1", buffer=5)

    # Adjust rotation for better view of the overall structure
    cmd.turn("y", 45)
    cmd.turn("x", -20)

def render_and_save():
    """Ray trace and save the image with enhanced shadows"""
    # Ray trace mode 1 = shadows enabled
    cmd.set("ray_trace_mode", 1)

    # Lighting settings for nice rendering with visible shadows
    cmd.set("ambient", 0.3)          # Lower ambient = darker shadows
    cmd.set("direct", 0.7)           # Higher direct light = more contrast
    cmd.set("specular", 0.4)         # Specular highlights
    cmd.set("shininess", 40)         # Shininess of highlights

    # Shadow settings - enhanced for visibility
    cmd.set("ray_shadows", 1)        # Enable shadows
    cmd.set("ray_shadow_decay_factor", 0.3)  # Shadow softness (higher = softer)
    cmd.set("ray_shadow_decay_range", 2.0)   # Shadow range

    # Light position for good shadow casting
    cmd.set("light_count", 2)
    cmd.set("light", "[-0.4, -0.4, -1.0]")  # Main light direction

    # Anti-aliasing
    cmd.set("antialias", 2)

    # Output size: 2 inch x 2 inch at 600 DPI = 1200 x 1200 pixels
    width_px = 1200
    height_px = 1200
    output_dpi = 600

    # Save session
    session_path = os.path.join(OUTPUT_DIR, "Figure6F_HYOU1_session.pse")
    cmd.save(session_path)
    print(f"Session saved: {session_path}")

    # Render with labels
    cmd.ray(width_px, height_px)
    png_path = os.path.join(OUTPUT_DIR, "Figure6F_HYOU1_structure.png")
    cmd.png(png_path, dpi=output_dpi)
    print(f"PNG saved: {png_path}")

    # Convert to PDF
    pdf_path = os.path.join(OUTPUT_DIR, "Figure6F_HYOU1_structure.pdf")
    try:
        import subprocess
        result = subprocess.run(
            ["sips", "-s", "format", "pdf", png_path, "--out", pdf_path],
            capture_output=True, text=True
        )
        if result.returncode == 0:
            print(f"PDF saved: {pdf_path}")
    except Exception as e:
        print(f"PDF conversion error: {e}")

def print_info():
    """Print site information"""
    print("\n" + "="*60)
    print("Figure 6F - HYOU1 (Q9Y4L1) - STRUCTURED REGION EXAMPLE")
    print("="*60)
    print(f"\nO-GlcNAc Site: T{SITE['residue']}")
    print(f"  - logFC: +{SITE['logFC']:.2f} (UP-regulated with Tunicamycin)")
    print(f"  - pLDDT: {SITE['pLDDT']} (High confidence = Structured)")
    print(f"  - Region: {SITE['region']}")
    print("\nVisualization:")
    print("  - Cartoon: Secondary structure colored by pLDDT")
    print("  - Surface: Transparent overlay (70% transparency)")
    print("  - Orange sphere: O-GlcNAc site T871")
    print("\npLDDT color scheme:")
    print("  - Dark blue (>90): Very high confidence")
    print("  - Light blue (70-90): High confidence")
    print("  - Yellow (50-70): Low confidence")
    print("  - Orange (<50): Very low confidence (IDR)")
    print("="*60 + "\n")

# ============================================
# Execute
# ============================================

if __name__ == "pymol" or __name__ == "__main__":
    print("Generating Figure 6F - HYOU1...")

    setup_pymol()
    load_structure()
    color_by_plddt()
    show_cartoon_and_surface()
    highlight_oglcnac_site()
    setup_view()
    print_info()
    render_and_save()

    print("\nDone! Adjust view manually if needed.")
