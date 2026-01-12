#!/bin/bash
# Setup script for StructureMap Python environment
# Used by Figure6.R for extracting structural features from AlphaFold predictions

# Create conda environment
echo "Creating structuremap_env conda environment..."
conda create -n structuremap_env python=3.9 -y

# Activate environment
echo "Activating environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate structuremap_env

# Install required packages
echo "Installing required packages..."
pip install structuremap pandas numpy scipy biopython

# Create AlphaFold directories
ALPHAFOLD_BASE="/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/reference"
echo ""
echo "Creating AlphaFold directories..."
mkdir -p "${ALPHAFOLD_BASE}/alphafold_cif"
mkdir -p "${ALPHAFOLD_BASE}/alphafold_pae"

# Print summary
echo ""
echo "=== Setup Complete ==="
echo "Environment created: structuremap_env"
echo ""
echo "AlphaFold directories created:"
echo "  CIF: ${ALPHAFOLD_BASE}/alphafold_cif/"
echo "  PAE: ${ALPHAFOLD_BASE}/alphafold_pae/"
echo ""
echo "Figure6.R is pre-configured to use these paths."
echo "Just run Figure6.R after this setup completes."
