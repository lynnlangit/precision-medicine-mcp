"""Spatial context extraction from spatial transcriptomics data.

Extracts spatial neighborhoods, computes distances, and identifies
relevant spatial relationships for quantum fidelity analysis.
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Optional, Any
from dataclasses import dataclass
from scipy.spatial import KDTree
from scipy.spatial.distance import pdist, squareform

try:
    import anndata as ad
    ANNDATA_AVAILABLE = True
except ImportError:
    ANNDATA_AVAILABLE = False
    ad = None


@dataclass
class SpatialNeighborhood:
    """Container for spatial neighborhood information.

    Attributes:
        cell_idx: Index of central cell
        neighbor_indices: Indices of neighboring cells
        distances: Distances to each neighbor
        cell_types: Cell type labels for neighbors
        coordinates: Spatial coordinates for neighbors
    """
    cell_idx: int
    neighbor_indices: List[int]
    distances: List[float]
    cell_types: List[str]
    coordinates: np.ndarray


class SpatialContextGenerator:
    """Extract spatial context from AnnData objects.

    Works with spatial transcriptomics data to identify:
    - K-nearest neighbors for each cell
    - Cells within a distance radius
    - Spatial patterns (clusters, boundaries, gradients)
    """

    def __init__(
        self,
        k_neighbors: int = 10,
        radius: Optional[float] = None,
        coordinate_keys: Tuple[str, str] = ("spatial_x", "spatial_y")
    ):
        """Initialize spatial context generator.

        Args:
            k_neighbors: Number of nearest neighbors to extract
            radius: Distance radius for neighborhood (if None, uses k_neighbors)
            coordinate_keys: Keys for spatial coordinates in adata.obs
        """
        self.k_neighbors = k_neighbors
        self.radius = radius
        self.coordinate_keys = coordinate_keys
        self.kdtree: Optional[KDTree] = None
        self.coordinates: Optional[np.ndarray] = None

    def fit(self, adata: Any) -> None:
        """Fit spatial index on AnnData coordinates.

        Args:
            adata: AnnData object with spatial coordinates
        """
        if not ANNDATA_AVAILABLE:
            raise ImportError("anndata not installed. Install with: pip install anndata")

        # Extract coordinates
        x_key, y_key = self.coordinate_keys

        if x_key not in adata.obs.columns or y_key not in adata.obs.columns:
            raise ValueError(
                f"Coordinate keys {self.coordinate_keys} not found in adata.obs. "
                f"Available keys: {list(adata.obs.columns)}"
            )

        self.coordinates = np.column_stack([
            adata.obs[x_key].values,
            adata.obs[y_key].values
        ])

        # Build KD-tree for efficient neighbor search
        self.kdtree = KDTree(self.coordinates)

        print(f"Fitted spatial index on {len(self.coordinates)} cells")

    def get_neighborhood(
        self,
        cell_idx: int,
        cell_type_labels: List[str],
        use_radius: bool = False
    ) -> SpatialNeighborhood:
        """Get spatial neighborhood for a cell.

        Args:
            cell_idx: Index of cell to get neighborhood for
            cell_type_labels: Cell type labels for all cells
            use_radius: If True, use radius-based neighbors; else use k-neighbors

        Returns:
            SpatialNeighborhood object
        """
        if self.kdtree is None or self.coordinates is None:
            raise RuntimeError("Must call fit() before getting neighborhoods")

        query_point = self.coordinates[cell_idx]

        if use_radius and self.radius is not None:
            # Radius-based neighbors
            indices = self.kdtree.query_ball_point(query_point, self.radius)
            # Remove self
            indices = [i for i in indices if i != cell_idx]
            # Compute distances
            if indices:
                distances = [
                    np.linalg.norm(self.coordinates[i] - query_point)
                    for i in indices
                ]
            else:
                distances = []
        else:
            # K-nearest neighbors
            k = min(self.k_neighbors + 1, len(self.coordinates))  # +1 for self
            distances, indices = self.kdtree.query(query_point, k=k)
            # Remove self (first entry)
            mask = indices != cell_idx
            indices = indices[mask].tolist()
            distances = distances[mask].tolist()

        # Get cell types for neighbors
        neighbor_cell_types = [cell_type_labels[i] for i in indices]

        # Get coordinates
        neighbor_coords = self.coordinates[indices]

        return SpatialNeighborhood(
            cell_idx=cell_idx,
            neighbor_indices=indices,
            distances=distances,
            cell_types=neighbor_cell_types,
            coordinates=neighbor_coords
        )

    def get_all_neighborhoods(
        self,
        cell_type_labels: List[str],
        cell_indices: Optional[List[int]] = None
    ) -> List[SpatialNeighborhood]:
        """Get neighborhoods for multiple cells.

        Args:
            cell_type_labels: Cell type labels for all cells
            cell_indices: Indices to get neighborhoods for (if None, use all)

        Returns:
            List of SpatialNeighborhood objects
        """
        if self.coordinates is None:
            raise RuntimeError("Must call fit() first")

        if cell_indices is None:
            cell_indices = list(range(len(self.coordinates)))

        neighborhoods = []
        for idx in cell_indices:
            neighborhood = self.get_neighborhood(idx, cell_type_labels)
            neighborhoods.append(neighborhood)

        return neighborhoods

    def compute_spatial_graph(
        self,
        cell_type_labels: List[str],
        include_edge_types: bool = True
    ) -> Dict[str, Any]:
        """Compute spatial proximity graph.

        Args:
            cell_type_labels: Cell type labels for all cells
            include_edge_types: If True, annotate edges with cell type pairs

        Returns:
            Dictionary with graph structure
        """
        if self.coordinates is None:
            raise RuntimeError("Must call fit() first")

        n_cells = len(self.coordinates)
        edges = []
        edge_distances = []
        edge_types = [] if include_edge_types else None

        # Get neighborhoods for all cells
        for cell_idx in range(n_cells):
            neighborhood = self.get_neighborhood(cell_idx, cell_type_labels)

            for neighbor_idx, distance in zip(
                neighborhood.neighbor_indices,
                neighborhood.distances
            ):
                # Add edge (undirected, so only add if cell_idx < neighbor_idx)
                if cell_idx < neighbor_idx:
                    edges.append((cell_idx, neighbor_idx))
                    edge_distances.append(distance)

                    if include_edge_types:
                        edge_types.append((
                            cell_type_labels[cell_idx],
                            cell_type_labels[neighbor_idx]
                        ))

        graph = {
            "n_cells": n_cells,
            "n_edges": len(edges),
            "edges": edges,
            "edge_distances": edge_distances,
            "edge_types": edge_types,
            "cell_type_labels": cell_type_labels
        }

        return graph

    def identify_tls_candidates(
        self,
        cell_type_labels: List[str],
        tls_marker_types: List[str],
        min_cluster_size: int = 20,
        max_distance: float = 100.0
    ) -> List[Dict[str, Any]]:
        """Identify potential tertiary lymphoid structures (TLS).

        TLS are characterized by dense clusters of immune cells
        (B cells, T cells, dendritic cells) in tumor tissue.

        Args:
            cell_type_labels: Cell type labels for all cells
            tls_marker_types: Cell types indicating TLS (e.g., ["B_cell", "T_cell"])
            min_cluster_size: Minimum cells to constitute TLS
            max_distance: Maximum distance for cells to be in same TLS

        Returns:
            List of TLS candidate dictionaries with cell indices and metadata
        """
        if self.coordinates is None:
            raise RuntimeError("Must call fit() first")

        # Find cells matching TLS marker types
        tls_cell_indices = [
            i for i, label in enumerate(cell_type_labels)
            if label in tls_marker_types
        ]

        if len(tls_cell_indices) < min_cluster_size:
            return []  # Not enough TLS marker cells

        # Get coordinates for TLS marker cells
        tls_coords = self.coordinates[tls_cell_indices]

        # Compute pairwise distances
        distances = squareform(pdist(tls_coords))

        # Find dense clusters using distance threshold
        # Simple approach: cells within max_distance form a cluster
        visited = set()
        clusters = []

        def dfs_cluster(start_idx, cluster):
            """DFS to find connected components."""
            if start_idx in visited:
                return
            visited.add(start_idx)
            cluster.append(tls_cell_indices[start_idx])

            # Find neighbors within max_distance
            neighbors = np.where(distances[start_idx] <= max_distance)[0]
            for neighbor_idx in neighbors:
                if neighbor_idx not in visited:
                    dfs_cluster(neighbor_idx, cluster)

        for i in range(len(tls_cell_indices)):
            if i not in visited:
                cluster = []
                dfs_cluster(i, cluster)
                if len(cluster) >= min_cluster_size:
                    clusters.append(cluster)

        # Create TLS candidate objects
        tls_candidates = []
        for cluster_idx, cluster_cells in enumerate(clusters):
            # Get cell type composition
            cluster_cell_types = [cell_type_labels[i] for i in cluster_cells]
            type_counts = pd.Series(cluster_cell_types).value_counts().to_dict()

            # Compute centroid
            cluster_coords = self.coordinates[cluster_cells]
            centroid = np.mean(cluster_coords, axis=0)

            # Compute spatial extent (max distance from centroid)
            centroid_distances = np.linalg.norm(
                cluster_coords - centroid,
                axis=1
            )
            spatial_extent = float(np.max(centroid_distances))

            tls_candidates.append({
                "tls_id": cluster_idx,
                "cell_indices": cluster_cells,
                "n_cells": len(cluster_cells),
                "cell_type_composition": type_counts,
                "centroid": centroid.tolist(),
                "spatial_extent": spatial_extent,
                "density": len(cluster_cells) / (np.pi * spatial_extent ** 2)
            })

        return tls_candidates

    def extract_features_for_adata(
        self,
        adata: Any,
        cell_type_key: str = "cell_type",
        gene_list: Optional[List[str]] = None,
        normalize: bool = True
    ) -> Tuple[np.ndarray, List[str]]:
        """Extract gene expression features from AnnData.

        Args:
            adata: AnnData object
            cell_type_key: Key in adata.obs for cell type labels
            gene_list: List of genes to extract (if None, use all)
            normalize: Whether to normalize features

        Returns:
            Tuple of (features_array, cell_type_labels)
            - features_array: (n_cells, n_genes)
            - cell_type_labels: List of cell type labels
        """
        if not ANNDATA_AVAILABLE:
            raise ImportError("anndata not installed")

        # Get cell type labels
        if cell_type_key not in adata.obs.columns:
            raise ValueError(f"Cell type key '{cell_type_key}' not found in adata.obs")

        cell_type_labels = adata.obs[cell_type_key].tolist()

        # Extract expression matrix
        if gene_list is not None:
            # Filter to specific genes
            gene_indices = [i for i, g in enumerate(adata.var_names) if g in gene_list]
            X = adata.X[:, gene_indices]
        else:
            X = adata.X

        # Convert to dense if sparse
        if hasattr(X, 'toarray'):
            X = X.toarray()

        # Normalize if requested
        if normalize:
            # Simple row normalization (per-cell)
            row_sums = X.sum(axis=1, keepdims=True)
            row_sums[row_sums == 0] = 1  # Avoid division by zero
            X = X / row_sums

        return X, cell_type_labels

    def get_summary(self) -> Dict[str, Any]:
        """Get summary of spatial context.

        Returns:
            Dictionary with summary information
        """
        if self.coordinates is None:
            return {"fitted": False}

        return {
            "fitted": True,
            "n_cells": len(self.coordinates),
            "k_neighbors": self.k_neighbors,
            "radius": self.radius,
            "coordinate_keys": self.coordinate_keys,
            "coordinate_range": {
                "x_min": float(np.min(self.coordinates[:, 0])),
                "x_max": float(np.max(self.coordinates[:, 0])),
                "y_min": float(np.min(self.coordinates[:, 1])),
                "y_max": float(np.max(self.coordinates[:, 1]))
            }
        }
