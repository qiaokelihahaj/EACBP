"""
EACBP Spatial Single-Cell Analytics Plane.
Provides native capabilities for spatial domain identification, spatial autocorrelation
(Moran's I / Geary's C), spatial differential expression, and ligand-receptor communication.
"""

from typing import Optional, List, Dict, Any
import numpy as np
import pandas as pd

from eacbp.capabilities.sc_data import SCData
from eacbp.capabilities.spatial.domain import (
    SpatialDomainCapability,
    validate_spatial_coordinates,
    build_spatial_neighborhood_graph,
    compute_spatially_smoothed_embedding,
    calculate_silhouette,
)
from eacbp.capabilities.spatial.autocorrelation import (
    SpatialDEGCapability,
    calculate_morans_i,
    calculate_gearys_c,
    benjamini_hochberg,
)
from eacbp.capabilities.spatial.cci import (
    CellCellCommunicationCapability,
    compute_spatial_cci,
    calculate_spatial_contact_density,
    CURATED_LIGAND_RECEPTOR_PAIRS,
)


def create_synthetic_spatial_ad_study(
    n_cells: int = 600,
    n_genes: int = 100,
    n_plaques: int = 3,
    random_seed: int = 42,
    grid_size: float = 1000.0,
) -> SCData:
    """
    Generates a realistic synthetic spatial transcriptomics dataset for AD research.
    
    Includes:
    - 2D spatial coordinates in .obsm['spatial']
    - Plaque foci generating localized microenvironments
    - Microglia, Astrocytes, Neurons, Oligodendrocytes
    - Known spatially localized markers (Apoe, Trem2, Clec7a)
    - Dispersed checkerboard markers
    - Random uniform background genes
    - Ligand-receptor communication pairs (Apoe-Trem2, App-Cd74, Cx3cl1-Cx3cr1)
    """
    np.random.seed(random_seed)

    # 1. Generate 2D spatial coordinates
    coords = np.random.uniform(0, grid_size, size=(n_cells, 2)).astype(np.float32)

    # Plaque centers
    plaque_centers = np.random.uniform(grid_size * 0.2, grid_size * 0.8, size=(n_plaques, 2))
    
    # Distance to closest plaque
    dists_to_plaques = np.min(
        np.linalg.norm(coords[:, None, :] - plaque_centers[None, :, :], axis=2),
        axis=1
    )
    is_plaque_adjacent = dists_to_plaques < (grid_size * 0.20)

    # 2. Cell types with spatial distribution
    cell_types = []
    for i in range(n_cells):
        if is_plaque_adjacent[i]:
            # Enriched for microglia and reactive astrocytes around plaques
            ct = np.random.choice(["Microglia", "Astrocytes", "Neurons", "Oligodendrocytes"], p=[0.50, 0.30, 0.15, 0.05])
        else:
            # Normal parenchyma
            ct = np.random.choice(["Microglia", "Astrocytes", "Neurons", "Oligodendrocytes"], p=[0.20, 0.30, 0.30, 0.20])
        cell_types.append(ct)
    cell_types = np.array(cell_types)

    # Conditions & mice
    conditions = ["AD" if is_plaque_adjacent[i] or np.random.rand() > 0.5 else "control" for i in range(n_cells)]
    mouse_ids = [f"AD_mouse_{(i % 4) + 1:02d}" if conditions[i] == "AD" else f"Ctrl_mouse_{(i % 4) + 1:02d}" for i in range(n_cells)]

    # 3. Gene names
    gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]
    marker_map = {
        0: "Apoe",       # High in plaque-adjacent microglia/astrocytes (Strong spatial positive autocorrelation)
        1: "Trem2",      # High in plaque-adjacent microglia (Strong spatial positive autocorrelation)
        2: "Clec7a",     # High in plaque-adjacent microglia (Strong spatial positive autocorrelation)
        3: "Cx3cr1",     # Homeostatic microglia marker
        4: "P2ry12",     # Homeostatic microglia marker (reduced near plaques)
        5: "Gfap",       # Astrocytes
        6: "Rbfox3",     # Neurons
        7: "Mog",        # Oligodendrocytes
        8: "App",        # Neuronal amyloid precursor protein (ligand for Cd74)
        9: "Cd74",       # Microglial receptor for App
        10: "Cx3cl1",    # Neuronal ligand
        11: "C3",        # Astrocytic / Microglial complement ligand
        12: "C3ar1",     # Microglial complement receptor
        13: "Spp1",      # Reactive microglia ligand
        14: "Cd44",      # Astrocyte / Microglia receptor
        15: "Dispersed_1", # Checkerboard pattern (Negative spatial autocorrelation)
        16: "Dispersed_2", # Checkerboard pattern
    }
    for idx, g in marker_map.items():
        if idx < n_genes:
            gene_names[idx] = g

    # 4. Base expression matrix
    base_rate = np.random.gamma(2.0, 1.0, size=(n_cells, n_genes))

    # Spatial proximity factor (higher near plaques)
    proximity_factor = np.exp(-dists_to_plaques / (grid_size * 0.18))

    for i in range(n_cells):
        ct = cell_types[i]
        prox = proximity_factor[i]
        
        # Base markers per cell type
        if ct == "Microglia":
            base_rate[i, 3] *= 4.0   # Cx3cr1
            base_rate[i, 4] *= 3.5   # P2ry12
            base_rate[i, 9] *= 3.0   # Cd74
            base_rate[i, 12] *= 3.0  # C3ar1
            base_rate[i, 4] *= max(0.1, 1.0 - 0.8 * prox) # P2ry12 downregulated
            
        elif ct == "Astrocytes":
            base_rate[i, 5] *= 6.0   # Gfap
            base_rate[i, 11] *= 4.0  # C3
            base_rate[i, 14] *= 3.5  # Cd44
            
        elif ct == "Neurons":
            base_rate[i, 6] *= 6.0   # Rbfox3
            base_rate[i, 8] *= 5.0   # App
            base_rate[i, 10] *= 5.0  # Cx3cl1
            
        elif ct == "Oligodendrocytes":
            base_rate[i, 7] *= 6.0   # Mog

        # Injected negative spatial autocorrelation (checkerboard frequency)
        x_bin = int(coords[i, 0] / (grid_size / 20.0)) % 2
        y_bin = int(coords[i, 1] / (grid_size / 20.0)) % 2
        if (x_bin + y_bin) % 2 == 0:
            base_rate[i, 15] *= 6.0
            base_rate[i, 16] *= 5.0
        else:
            base_rate[i, 15] *= 0.2
            base_rate[i, 16] *= 0.2

    # Add strong plaque-associated spatial hotspots
    base_rate[:, 0] += 20.0 * proximity_factor   # Apoe
    base_rate[:, 1] += 16.0 * proximity_factor   # Trem2
    base_rate[:, 2] += 16.0 * proximity_factor   # Clec7a
    base_rate[:, 13] += 14.0 * proximity_factor  # Spp1

    # Sample Poisson counts
    X = np.random.poisson(base_rate).astype(np.float32)

    obs = pd.DataFrame({
        "cell_id": [f"cell_{i:05d}" for i in range(n_cells)],
        "mouse_id": mouse_ids,
        "condition": conditions,
        "cell_type": cell_types,
        "cell_type_ground_truth": cell_types,
        "x_coord": coords[:, 0],
        "y_coord": coords[:, 1],
        "dist_to_nearest_plaque": dists_to_plaques,
        "is_plaque_adjacent": is_plaque_adjacent,
        "n_counts": X.sum(axis=1),
        "n_genes": (X > 0).sum(axis=1),
    })

    var = pd.DataFrame({
        "gene_name": gene_names,
        "n_cells": (X > 0).sum(axis=0),
    })

    obsm = {
        "spatial": coords,
        "X_pca": base_rate[:, :min(20, n_genes)].astype(np.float32),
    }

    return SCData(X=X, obs=obs, var=var, obsm=obsm)


__all__ = [
    "SpatialDomainCapability",
    "SpatialDEGCapability",
    "CellCellCommunicationCapability",
    "validate_spatial_coordinates",
    "build_spatial_neighborhood_graph",
    "compute_spatially_smoothed_embedding",
    "calculate_silhouette",
    "calculate_morans_i",
    "calculate_gearys_c",
    "benjamini_hochberg",
    "compute_spatial_cci",
    "calculate_spatial_contact_density",
    "create_synthetic_spatial_ad_study",
    "CURATED_LIGAND_RECEPTOR_PAIRS",
]
