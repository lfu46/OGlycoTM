"""
UMAP Analysis for O-GlcNAc Proteomics Data
Called from R via reticulate
"""

import numpy as np
import pandas as pd
from umap import UMAP
from sklearn.preprocessing import StandardScaler


def run_umap(intensity_matrix, n_neighbors=5, min_dist=0.3, n_components=2,
             metric='euclidean', random_state=42):
    """
    Run UMAP on intensity matrix.

    Parameters:
    -----------
    intensity_matrix : numpy array or pandas DataFrame
        Samples as rows, proteins/features as columns
    n_neighbors : int
        Number of neighbors for UMAP (default 5, good for small datasets)
    min_dist : float
        Minimum distance between points in UMAP space (default 0.3)
    n_components : int
        Number of UMAP dimensions (default 2)
    metric : str
        Distance metric (default 'euclidean')
    random_state : int
        Random seed for reproducibility

    Returns:
    --------
    numpy array : UMAP coordinates (n_samples x n_components)
    """
    # Convert to numpy if needed
    if isinstance(intensity_matrix, pd.DataFrame):
        intensity_matrix = intensity_matrix.values

    # Standardize the data (important for UMAP)
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(intensity_matrix)

    # Run UMAP
    umap_model = UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        n_components=n_components,
        metric=metric,
        random_state=random_state
    )

    umap_coords = umap_model.fit_transform(scaled_data)

    return umap_coords


def run_umap_with_metadata(intensity_df, sample_names, cell_types, conditions,
                           n_neighbors=5, min_dist=0.3, random_state=42):
    """
    Run UMAP and return a DataFrame with coordinates and metadata.

    Parameters:
    -----------
    intensity_df : pandas DataFrame
        Samples as rows, proteins as columns
    sample_names : list
        Sample identifiers
    cell_types : list
        Cell type for each sample
    conditions : list
        Condition (Tuni/Ctrl) for each sample
    n_neighbors : int
        UMAP parameter
    min_dist : float
        UMAP parameter
    random_state : int
        Random seed

    Returns:
    --------
    pandas DataFrame with UMAP1, UMAP2, Sample, CellType, Condition columns
    """
    # Run UMAP
    umap_coords = run_umap(
        intensity_df,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        random_state=random_state
    )

    # Create result DataFrame
    result_df = pd.DataFrame({
        'UMAP1': umap_coords[:, 0],
        'UMAP2': umap_coords[:, 1],
        'Sample': sample_names,
        'CellType': cell_types,
        'Condition': conditions
    })

    return result_df


def run_umap_proteins(logfc_matrix, protein_ids, n_neighbors=15, min_dist=0.1,
                       random_state=42):
    """
    Run UMAP on protein logFC matrix for protein-centric visualization.

    Parameters:
    -----------
    logfc_matrix : pandas DataFrame or numpy array
        Proteins as rows, cell types (logFC values) as columns
    protein_ids : list
        Protein identifiers
    n_neighbors : int
        UMAP parameter (default 15, good for larger datasets)
    min_dist : float
        UMAP parameter (default 0.1)
    random_state : int
        Random seed

    Returns:
    --------
    pandas DataFrame with UMAP1, UMAP2, Protein.ID columns
    """
    # Run UMAP
    umap_coords = run_umap(
        logfc_matrix,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        random_state=random_state
    )

    # Create result DataFrame
    result_df = pd.DataFrame({
        'UMAP1': umap_coords[:, 0],
        'UMAP2': umap_coords[:, 1],
        'Protein.ID': protein_ids
    })

    return result_df


if __name__ == "__main__":
    # Test with dummy data
    np.random.seed(42)

    # Create dummy data: 18 samples (3 cell types x 2 conditions x 3 replicates)
    n_proteins = 100
    n_samples = 18

    dummy_data = np.random.randn(n_samples, n_proteins)

    sample_names = []
    cell_types = []
    conditions = []

    for cell in ['HEK293T', 'HepG2', 'Jurkat']:
        for cond in ['Tuni', 'Ctrl']:
            for rep in [1, 2, 3]:
                sample_names.append(f"{cell}_{cond}_{rep}")
                cell_types.append(cell)
                conditions.append(cond)

    result = run_umap_with_metadata(
        pd.DataFrame(dummy_data),
        sample_names,
        cell_types,
        conditions
    )

    print("UMAP Result:")
    print(result)
