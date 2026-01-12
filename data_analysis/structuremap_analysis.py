"""
StructureMap Analysis for O-GlcNAc Site Feature Extraction

This script uses the StructureMap Python package to extract structural features
from AlphaFold predictions for O-GlcNAc sites.

Features extracted:
- pPSE_24 (nAA_24_180_pae): Full sphere exposure, 24A radius, 180 angle - for IDR prediction
- pPSE_12 (nAA_12_70_pae): Part sphere exposure, 12A radius, 70 angle - for side chain accessibility
- Secondary structure: HELX, STRN, BEND, TURN, or unstructured
- IDR status: Intrinsically disordered region (smoothed pPSE_24 <= 34.27)

Usage:
    Run via reticulate from R, or directly in Python.
"""

import os
import warnings
import numpy as np
import pandas as pd

# Suppress FutureWarnings from structuremap/pandas
warnings.filterwarnings('ignore', category=FutureWarning)
from structuremap.processing import (
    download_alphafold_cif,
    download_alphafold_pae,
    format_alphafold_data,
    annotate_accessibility,
    get_smooth_score
)


def download_alphafold_files(protein_ids, cif_folder, pae_folder):
    """
    Download AlphaFold CIF and PAE files for a list of proteins.

    Parameters:
    -----------
    protein_ids : list
        List of UniProt accession IDs
    cif_folder : str
        Output folder for CIF files
    pae_folder : str
        Output folder for PAE files

    Returns:
    --------
    dict : Dictionary with download statistics
    """
    # Create output directories if they don't exist
    os.makedirs(cif_folder, exist_ok=True)
    os.makedirs(pae_folder, exist_ok=True)

    # Clean protein IDs (remove any whitespace)
    protein_ids = [p.strip() for p in protein_ids if isinstance(p, str) and len(p.strip()) > 0]

    print(f"Downloading AlphaFold files for {len(protein_ids)} proteins...")
    print(f"CIF folder: {cif_folder}")
    print(f"PAE folder: {pae_folder}")

    # Download CIF files
    print("\nDownloading CIF files...")
    cif_result = download_alphafold_cif(
        proteins=protein_ids,
        out_folder=cif_folder
    )

    # Download PAE files
    print("\nDownloading PAE files...")
    pae_result = download_alphafold_pae(
        proteins=protein_ids,
        out_folder=pae_folder
    )

    # Count downloaded files
    cif_files = [f for f in os.listdir(cif_folder) if f.endswith('.cif')]
    pae_files = [f for f in os.listdir(pae_folder) if f.endswith('.json')]

    stats = {
        'total_proteins': len(protein_ids),
        'cif_downloaded': len(cif_files),
        'pae_downloaded': len(pae_files)
    }

    print(f"\nDownload complete:")
    print(f"  CIF files: {stats['cif_downloaded']}/{stats['total_proteins']}")
    print(f"  PAE files: {stats['pae_downloaded']}/{stats['total_proteins']}")

    return stats


def extract_structural_features(protein_ids, cif_folder, pae_folder):
    """
    Extract structural features (pPSE and secondary structure) for proteins.

    Parameters:
    -----------
    protein_ids : list
        List of UniProt accession IDs
    cif_folder : str
        Folder containing CIF files
    pae_folder : str
        Folder containing PAE files

    Returns:
    --------
    pandas.DataFrame : DataFrame with structural features per residue
    """
    # Clean protein IDs
    protein_ids = [p.strip() for p in protein_ids if isinstance(p, str) and len(p.strip()) > 0]

    print(f"Extracting structural features for {len(protein_ids)} proteins...")

    # Step 1: Format AlphaFold data (extract coordinates and secondary structure)
    print("Step 1: Formatting AlphaFold data...")
    alphafold_df = format_alphafold_data(
        directory=cif_folder,
        protein_ids=protein_ids
    )

    if alphafold_df is None or len(alphafold_df) == 0:
        print("Warning: No AlphaFold data could be formatted.")
        return None

    print(f"  Formatted {len(alphafold_df)} residues from {alphafold_df['protein_id'].nunique()} proteins")

    # Step 2: Annotate full sphere exposure (pPSE_24 for IDR prediction)
    # Parameters: 24 Angstrom radius, 180 degree angle
    # Exactly matching tutorial cell-24
    print("Step 2: Calculating pPSE_24 (full sphere, 24A, 180deg)...")
    full_sphere_exposure = annotate_accessibility(
        df=alphafold_df,
        max_dist=24,
        max_angle=180,
        error_dir=pae_folder
    )

    # Step 3: Annotate part sphere exposure (pPSE_12 for side chain accessibility)
    # Parameters: 12 Angstrom radius, 70 degree angle
    # Exactly matching tutorial cell-29
    print("Step 3: Calculating pPSE_12 (part sphere, 12A, 70deg)...")
    half_sphere_exposure = annotate_accessibility(
        df=alphafold_df,
        max_dist=12,
        max_angle=70,
        error_dir=pae_folder
    )

    # Combine results - exactly matching tutorial cell-27 and cell-32
    print("Step 4: Combining results...")

    # Merge full sphere exposure with alphafold_annotation
    alphafold_accessibility = alphafold_df.merge(
        full_sphere_exposure, how='left', on=['protein_id', 'AA', 'position']
    )

    # Merge half sphere exposure
    alphafold_accessibility = alphafold_accessibility.merge(
        half_sphere_exposure, how='left', on=['protein_id', 'AA', 'position']
    )

    # Step 5: Calculate smoothed pPSE_24 and IDR status
    # Exactly matching tutorial cell-35
    print("Step 5: Calculating smoothed pPSE_24 for IDR prediction...")

    alphafold_accessibility_smooth = get_smooth_score(
        alphafold_accessibility,
        np.array(['nAA_24_180_pae']),
        [10]
    )

    # Predict IDR based on threshold (smoothed pPSE_24 <= 34.27)
    # Exactly matching tutorial cell-37
    alphafold_accessibility_smooth['is_IDR'] = np.where(
        alphafold_accessibility_smooth['nAA_24_180_pae_smooth10'] <= 34.27, 1, 0
    )

    # Use this as result_df
    result_df = alphafold_accessibility_smooth

    # Rename columns for clarity at the end
    result_df = result_df.rename(columns={
        'nAA_24_180_pae': 'pPSE_24',
        'nAA_12_70_pae': 'pPSE_12',
        'nAA_24_180_pae_smooth10': 'pPSE_24_smooth10',
        'quality': 'pLDDT'
    })

    # Step 6: Map secondary structure to simplified categories
    print("Step 6: Mapping secondary structure...")

    # Secondary structure mapping
    # From CIF file: HELX_P (alpha helix), STRN (strand), BEND, TURN
    # Tutorial shows columns: BEND, HELX, STRN, TURN, unstructured (cell-21)
    ss_mapping = {
        'HELX_P': 'HELX',
        'STRN': 'STRN',
        'BEND': 'BEND',
        'TURN': 'TURN'
    }

    if 'secondary_structure' in result_df.columns:
        result_df['secondary_structure_simple'] = result_df['secondary_structure'].map(
            lambda x: ss_mapping.get(x, 'unstructured') if pd.notna(x) else 'unstructured'
        )
    else:
        result_df['secondary_structure_simple'] = 'unstructured'

    # Select final columns
    final_columns = [
        'protein_id',
        'AA',  # Amino acid
        'position',  # Position in protein (1-indexed)
        'pLDDT',  # AlphaFold confidence (renamed from 'quality')
        'secondary_structure',
        'secondary_structure_simple',
        'pPSE_24',
        'pPSE_12',
        'pPSE_24_smooth10',
        'is_IDR'
    ]

    # Filter to only available columns
    available_columns = [col for col in final_columns if col in result_df.columns]
    result_df = result_df[available_columns].reset_index(drop=True)

    print(f"\nExtraction complete: {len(result_df)} residues with structural features")
    print(f"  Columns: {list(result_df.columns)}")

    return result_df


def get_site_structural_features(structural_df, sites_df):
    """
    Extract structural features for specific O-GlcNAc sites.

    Parameters:
    -----------
    structural_df : pandas.DataFrame
        DataFrame with structural features (from extract_structural_features)
    sites_df : pandas.DataFrame
        DataFrame with O-GlcNAc sites (must have 'Protein.ID' and 'site_number' columns)

    Returns:
    --------
    pandas.DataFrame : Input sites_df with additional structural feature columns
    """
    # Ensure column names are correct
    if 'Protein.ID' in sites_df.columns:
        sites_df = sites_df.rename(columns={'Protein.ID': 'protein_id'})
    if 'site_number' in sites_df.columns:
        sites_df = sites_df.rename(columns={'site_number': 'position'})

    # Convert position to integer for merging
    sites_df['position'] = sites_df['position'].astype(int)
    structural_df['position'] = structural_df['position'].astype(int)

    # Merge structural features with sites
    result = sites_df.merge(
        structural_df,
        on=['protein_id', 'position'],
        how='left'
    )

    # Report coverage
    total_sites = len(sites_df)
    sites_with_features = result['pPSE_24'].notna().sum()

    print(f"\nStructural feature coverage: {sites_with_features}/{total_sites} sites ({100*sites_with_features/total_sites:.1f}%)")

    return result


def main():
    """
    Main function for standalone execution (testing).
    """
    # Example usage - to be run from R via reticulate
    print("StructureMap Analysis Module")
    print("Use via reticulate in R or import functions directly.")
    print("\nAvailable functions:")
    print("  - download_alphafold_files(protein_ids, cif_folder, pae_folder)")
    print("  - extract_structural_features(protein_ids, cif_folder, pae_folder)")
    print("  - get_site_structural_features(structural_df, sites_df)")


if __name__ == "__main__":
    main()
