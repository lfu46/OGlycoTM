from pymol import cmd, stored
import os

STRUCTURE_PATH = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/alphafold_structures/AF-P02786-F1-model_v6.pdb"
DEFAULT_OUTPUT_DIR = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/Figure6_OGalNAc"
FALLBACK_OUTPUT_DIR = "/Users/longpingfu/Downloads/OGlycoTM/data_analysis/Figure6_OGalNAc_preview/secretory_downregulated"


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

def setup_pymol():
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
    cmd.load(STRUCTURE_PATH, "TFRC")
    cmd.hide("everything")

def color_by_plddt():
    cmd.set_color("plddt_very_high", [0/255, 83/255, 214/255])
    cmd.set_color("plddt_high", [101/255, 203/255, 243/255])
    cmd.set_color("plddt_low", [255/255, 219/255, 19/255])
    cmd.set_color("plddt_very_low", [255/255, 125/255, 69/255])
    stored.plddt_colors = {}
    cmd.iterate("TFRC and name CA", "stored.plddt_colors[resi] = b")
    for resi, plddt in stored.plddt_colors.items():
        if plddt > 90: cmd.color("plddt_very_high", f"TFRC and resi {resi}")
        elif plddt > 70: cmd.color("plddt_high", f"TFRC and resi {resi}")
        elif plddt > 50: cmd.color("plddt_low", f"TFRC and resi {resi}")
        else: cmd.color("plddt_very_low", f"TFRC and resi {resi}")

def show_cartoon_and_surface():
    cmd.show("cartoon", "TFRC")
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.show("surface", "TFRC")
    cmd.set("transparency", 0.7, "TFRC")
    cmd.set("surface_quality", 1)

def highlight_sites():
    # T104: O-GalNAc site, logFC=-3.17, protects against ectodomain shedding
    cmd.select("site_T104", "TFRC and resi 104")
    cmd.show("spheres", "site_T104")
    cmd.color("orange", "site_T104")  # orange = downregulated
    cmd.set("sphere_scale", 2.0, "site_T104")
    cmd.set("sphere_transparency", 0, "site_T104")

def setup_view():
    cmd.orient("TFRC")
    cmd.zoom("TFRC", buffer=5)
    cmd.turn("y", 45)
    cmd.turn("x", -20)

def render_and_save():
    cmd.set("ray_trace_mode", 1)
    cmd.set("ambient", 0.3)
    cmd.set("direct", 0.7)
    cmd.set("specular", 0.4)
    cmd.set("shininess", 40)
    cmd.set("ray_shadows", 1)
    cmd.set("ray_shadow_decay_factor", 0.3)
    cmd.set("ray_shadow_decay_range", 2.0)
    cmd.set("light_count", 2)
    cmd.set("light", "[-0.4, -0.4, -1.0]")
    cmd.set("antialias", 2)
    cmd.ray(1200, 1200)
    png_path = os.path.join(OUTPUT_DIR, "OGalNAc_TFRC_T104_structure.png")
    cmd.png(png_path, dpi=600)
    print(f"PNG saved: {png_path}")
    pdf_path = os.path.join(OUTPUT_DIR, "OGalNAc_TFRC_T104_structure.pdf")
    import subprocess
    subprocess.run(["sips", "-s", "format", "pdf", png_path, "--out", pdf_path], capture_output=True)
    print(f"PDF saved: {pdf_path}")
    session_path = os.path.join(OUTPUT_DIR, "OGalNAc_TFRC_T104_session.pse")
    cmd.save(session_path)

if __name__ == "pymol" or __name__ == "__main__":
    setup_pymol()
    load_structure()
    color_by_plddt()
    show_cartoon_and_surface()
    highlight_sites()
    setup_view()
    render_and_save()
    print("Done!")
