"""
Figure 6F - Compare Different Ray Trace Modes
Generates examples of each ray_trace_mode for comparison
"""

from pymol import cmd, stored
import os

# ============================================
# Configuration
# ============================================

STRUCTURE_PATH = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/alphafold_structures/AF-Q9NYJ8-F1-model_v6.pdb"
OUTPUT_DIR = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/Figure6"

# O-GlcNAc site data (HEK293T)
SITES = {
    29: {"logFC": 0.824, "region": "Structured", "pLDDT": 88.2},
    359: {"logFC": -0.272, "region": "IDR", "pLDDT": 31.8},
    460: {"logFC": -1.54, "region": "IDR", "pLDDT": 53.3}
}

COLOR_UP = "tv_orange"
COLOR_DOWN = "cyan"

# ============================================
# Setup Functions
# ============================================

def setup_pymol():
    cmd.bg_color("white")
    cmd.set("antialias", 2)
    cmd.set("ray_opaque_background", 1)
    cmd.set("depth_cue", 0)
    cmd.set("fog", 0)

def load_structure():
    cmd.load(STRUCTURE_PATH, "TAB2")
    cmd.hide("everything")
    cmd.show("cartoon", "TAB2")

def color_by_plddt():
    cmd.set_color("plddt_very_high", [0/255, 83/255, 214/255])
    cmd.set_color("plddt_high", [101/255, 203/255, 243/255])
    cmd.set_color("plddt_low", [255/255, 219/255, 19/255])
    cmd.set_color("plddt_very_low", [255/255, 125/255, 69/255])

    stored.plddt_colors = {}
    cmd.iterate("TAB2 and name CA", "stored.plddt_colors[resi] = b")

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
    for resi, data in SITES.items():
        site_name = f"site_S{resi}"
        cmd.select(site_name, f"TAB2 and resi {resi}")
        cmd.show("spheres", f"{site_name} and (name CB or name OG)")
        if data["logFC"] > 0:
            cmd.color(COLOR_UP, f"{site_name} and (name CB or name OG)")
        else:
            cmd.color(COLOR_DOWN, f"{site_name} and (name CB or name OG)")
        cmd.set("sphere_scale", 1.2, site_name)
    cmd.select("all_sites", " or ".join([f"site_S{r}" for r in SITES.keys()]))

def setup_view():
    cmd.orient("TAB2")
    cmd.turn("y", 45)
    cmd.turn("x", -15)
    cmd.zoom("TAB2", buffer=8)

# ============================================
# Ray Trace Mode Definitions
# ============================================

RAY_MODES = {
    "mode0_default": {
        "ray_trace_mode": 0,
        "description": "Mode 0: Default - Clean, simple rendering",
        "settings": {
            "ambient": 0.4,
            "direct": 0.5,
            "specular": 0.2,
            "ray_shadows": 0,
        }
    },
    "mode1_shadows": {
        "ray_trace_mode": 1,
        "description": "Mode 1: Shadows - Realistic with soft shadows",
        "settings": {
            "ambient": 0.35,
            "direct": 0.6,
            "specular": 0.3,
            "ray_shadows": 1,
            "ray_shadow_decay_factor": 0.1,
        }
    },
    "mode2_outline": {
        "ray_trace_mode": 2,
        "description": "Mode 2: Black Outline - Publication-friendly, clear edges",
        "settings": {
            "ambient": 0.5,
            "direct": 0.4,
            "specular": 0.1,
            "ray_shadows": 0,
            "ray_trace_disco_factor": 1.0,
        }
    },
    "mode3_quantized": {
        "ray_trace_mode": 3,
        "description": "Mode 3: Quantized - Posterized/cartoon look",
        "settings": {
            "ambient": 0.5,
            "direct": 0.5,
            "specular": 0.1,
            "ray_shadows": 0,
        }
    },
}

def render_mode(mode_name, mode_config):
    """Render structure with specific ray trace mode"""
    print(f"\nRendering {mode_config['description']}...")

    # Apply mode settings
    cmd.set("ray_trace_mode", mode_config["ray_trace_mode"])
    for setting, value in mode_config["settings"].items():
        cmd.set(setting, value)

    # Common quality settings
    cmd.set("antialias", 2)
    cmd.set("ray_trace_frames", 1)
    cmd.set("ray_opaque_background", 1)

    # Output size: 2 inch x 2 inch at 600 DPI
    width_px = 1200
    height_px = 1200

    # Ray trace and save
    cmd.ray(width_px, height_px)
    output_path = os.path.join(OUTPUT_DIR, f"Figure6F_TAB2_{mode_name}.png")
    cmd.png(output_path, dpi=600)
    print(f"Saved: {output_path}")

    return output_path

# ============================================
# Main
# ============================================

if __name__ == "pymol" or __name__ == "__main__":
    print("="*60)
    print("Generating Ray Trace Mode Comparison")
    print("="*60)

    # Setup
    setup_pymol()
    load_structure()
    color_by_plddt()
    highlight_oglcnac_sites()
    setup_view()

    # Hide labels for clean comparison
    cmd.hide("labels")

    # Render each mode
    for mode_name, mode_config in RAY_MODES.items():
        render_mode(mode_name, mode_config)

    print("\n" + "="*60)
    print("All modes rendered! Compare the following files:")
    print("="*60)
    for mode_name, mode_config in RAY_MODES.items():
        print(f"  {mode_name}.png - {mode_config['description']}")
    print("="*60)
