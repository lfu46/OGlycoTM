# Figure 6F - TAB2 (Q9NYJ8) 3D Structure with O-GlcNAc Sites
# PyMOL script for visualizing AlphaFold structure colored by pLDDT
# with O-GlcNAc sites highlighted as spheres

# Initialize settings
bg_color white
set ray_opaque_background, 1
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, 0

# Load AlphaFold structure
load /Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/alphafold_structures/AF-Q9NYJ8-F1-model_v6.pdb, TAB2

# Show as cartoon representation
hide everything
show cartoon, TAB2

# Color by pLDDT (stored in B-factor column)
# AlphaFold pLDDT color scheme: blue (high confidence) to red (low confidence)
spectrum b, blue_white_red, TAB2, minimum=50, maximum=90

# Alternative: Use AlphaFold official color scheme
# Very high (pLDDT > 90): dark blue
# High (70 < pLDDT < 90): light blue
# Low (50 < pLDDT < 70): yellow
# Very low (pLDDT < 50): orange/red

# Select O-GlcNAc sites
select site_S29, TAB2 and resi 29
select site_S359, TAB2 and resi 359
select site_S460, TAB2 and resi 460
select all_sites, site_S29 or site_S359 or site_S460

# Show sites as spheres (serine side chain)
show spheres, all_sites and (name CB or name OG)

# Color sites by fold change direction
# S29: logFC = 0.824 (UP - orange/salmon color)
# S359: logFC = -0.272 (DOWN - cyan)
# S460: logFC = -1.54 (DOWN - cyan)
color tv_orange, site_S29 and (name CB or name OG)
color cyan, site_S359 and (name CB or name OG)
color cyan, site_S460 and (name CB or name OG)

# Adjust sphere size
set sphere_scale, 1.5, all_sites

# Set up labels
# Note: Labels are added as pseudoatoms for better positioning
pseudoatom label_S29, pos=[30, 0, 0], label="S29\nFC=0.82"
pseudoatom label_S359, pos=[0, 0, 0], label="S359\nFC=-0.27"
pseudoatom label_S460, pos=[0, 0, 0], label="S460\nFC=-1.54"

# Get actual positions and offset labels
# These will be adjusted based on the actual structure orientation
iterate_state 1, site_S29 and name CA, stored.s29_pos = [x, y, z]
iterate_state 1, site_S359 and name CA, stored.s359_pos = [x, y, z]
iterate_state 1, site_S460 and name CA, stored.s460_pos = [x, y, z]

# Label styling
set label_font_id, 7
set label_size, 18
set label_color, black
set label_outline_color, white
set label_position, [2, 2, 2]

# Add labels directly to residues
label site_S29 and name CA, "S29\nFC=0.82"
label site_S359 and name CA, "S359\nFC=-0.27"
label site_S460 and name CA, "S460\nFC=-1.54"

# Clean up pseudoatoms (use direct labels instead)
delete label_S29
delete label_S359
delete label_S460

# Orient the molecule for best view
orient TAB2
turn y, 30
turn x, -20

# Zoom to fit
zoom TAB2, buffer=5

# Set up for high-quality output
set ray_trace_frames, 1
set ray_trace_gain, 0.1
set ambient, 0.4
set direct, 0.3
set specular, 0.15

# Save session for manual adjustment if needed
save /Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/Figure6F_TAB2_session.pse

# Ray trace and save image
ray 2400, 2400
png /Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/Figure6F_TAB2_structure.png, dpi=300

print "Figure 6F generated successfully!"
print "Output: /Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/Figure6F_TAB2_structure.png"
