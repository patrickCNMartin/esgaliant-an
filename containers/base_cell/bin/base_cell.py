import cellxgene_census
import numpy as np
import numpy.typing as npt
import tiledbsoma as soma
import zarr
import os
from typing import Literal

try:
    from scipy.stats import wasserstein_distance
except ImportError:
    wasserstein_distance = None


def compute_chunk_size(
    n_cells: int = 10000,
    n_genes: int = 2000,
    max_chunk_size: int = 10000,
    min_chunk_size: int = 1000,
):
    """
    Compute chunk sizes for zarr arrays.
    
    Parameters:
    -----------
    n_cells : int
        Number of cells (default: 10000)
    n_genes : int
        Number of genes (default: 2000)
    max_chunk_size : int
        Maximum chunk size (default: 10000)
    min_chunk_size : int
        Minimum chunk size (default: 1000)
    
    Returns:
    --------
    tuple[int, int]
        (gene_chunk_size, cell_chunk_size)
    """
    if n_genes < min_chunk_size:
        raise ValueError(
            "Number of genes should not be smaller than min chunk size."
        )
    if n_cells < min_chunk_size:
        raise ValueError(
            "Number of cells should not be smaller than min chunk size."
        )
    gene_chunk_size = min(max_chunk_size, n_genes)
    cell_chunk_size = min(max_chunk_size, n_cells)
    return gene_chunk_size, cell_chunk_size


def base_cell(
    atlas: str = "cellxgene",
    organism: str = "mus_musculus",
    cell_types: None | list[str] = None,
):
    if atlas == "cellxgene":
        gene_set = get_gene_set(atlas, organism)
        base_cell = base_cellxgene(gene_set, organism, cell_types)
        return gene_set, base_cell
    else:
        raise ValueError(f"Unknown atlas: {atlas}")


def get_gene_set(
    atlas: str = "cellxgene",
    organism: str = "mus_musculus",
):
    if atlas == "cellxgene":
        with cellxgene_census.open_soma(census_version="2025-01-30") as census:
            gene_metadata = cellxgene_census.get_var(
                census, 
                organism,
                column_names=["feature_name", "feature_id"])
            genes = set(gene_metadata['feature_name'].tolist())
        return genes
    else:
        raise ValueError(f"Unknown atlas: {atlas}")


def base_cellxgene(
    gene_set: list[str],
    organism: str = "mus_musculus",
    cell_types: list[str] | None = None,
) -> npt.NDArray[np.float64]:
    """
    Compute mean expression for a gene set, optionally filtered by cell types.

    Parameters:
    -----------
    gene_set : list[str]
        Compulsory. List of gene names to compute mean expression for.
    organism : str
        Organism name (default: "mus_musculus")
    cell_types : list[str] | None
        Optional. List of cell type values to filter on. If None, uses all cell types.

    Returns:
    --------
    npt.NDArray[np.float64]
        Mean expression values for the gene set across selected cells.
    """
    # connect to cellxgene server
    with cellxgene_census.open_soma(census_version="2025-01-30") as census:
        org = census["census_data"][organism]

        # Build obs_query with primary data filter
        obs_filter = "is_primary_data==True"

        # Add cell type filter if provided
        if cell_types is not None:
            cell_type_filter = " and ".join(
                [f'cell_type=="{ct}"' for ct in cell_types]
            )
            obs_filter = f"({obs_filter}) && ({cell_type_filter})"

        # Always filter primary tissue to remove duplicate cells
        with org.axis_query(
            measurement_name="RNA",
            obs_query=soma.AxisQuery(value_filter=obs_filter),
            var_query=soma.AxisQuery(
                value_filter=" or ".join(
                    [f'feature_name=="{gene}"' for gene in gene_set]
                )
            ),
        ) as query:
            # Get variable information
            var_df = query.var().concat().to_pandas()
            n_vars = len(var_df)
            n_obs = query.n_obs

            print(f"Filtering for {n_obs} cells and {n_vars} genes")

            # Initialize accumulators for per-gene means
            gene_sum = np.zeros((n_vars,), dtype=np.float64)
            gene_count = np.zeros((n_vars,), dtype=np.int64)

            # Get indexer to map soma_joinid to positional indices
            indexer = query.indexer

            # Stream through X data in batches
            for chunk_idx, arrow_tbl in enumerate(query.X("raw").tables()):
                print(f"Processing chunk {chunk_idx + 1}...")

                # Get positional indices for genes (var dimension)
                var_pos = indexer.by_var(arrow_tbl["soma_dim_1"])
                # Get the data values
                data = arrow_tbl["soma_data"].to_numpy()

                # Accumulate sums and counts per gene
                np.add.at(gene_sum, var_pos, data)
                np.add.at(gene_count, var_pos, 1)

                print(
                    f"  Chunk {chunk_idx + 1} complete: {len(data)} values processed"
                )

            # Compute final means
            gene_mean = np.divide(
                gene_sum,
                n_obs,
                where=(gene_count > 0),
                out=np.zeros_like(gene_sum),
            )
    return gene_mean


def cell_type_base(
    zarr_store,
    atlas: str = "cellxgene",
    organism: str = "mus_musculus",
    keep_all: bool = False,  # might remove this option
    max_chunk_size: int = 5000,
    min_chunk_size: int = 1000,
):
    """
    Compute mean expression per cell type and store in zarr store.

    Parameters:
    -----------
    zarr_store : str | zarr.Group
        Path to zarr store or zarr Group object
    atlas : str
        Atlas name (default: "cellxgene")
    organism : str
        Organism name (default: "mus_musculus")
    keep_all : bool
        If True, keep all cell types. If False, use unique cell types only.
    max_chunk_size : int
        Maximum chunk size for zarr arrays (default: 5000)
    min_chunk_size : int
        Minimum chunk size for zarr arrays (default: 1000)

    Returns:
    --------
    None
        Results are stored in the zarr store.
    """
    # Open zarr store if path is provided
    if isinstance(zarr_store, str):
        store = zarr.storage.LocalStore(zarr_store, read_only=False)
        root = zarr.group(store)
    else:
        root = zarr_store

    if atlas == "cellxgene":
        # Get the atlas group
        atlas_group_name = f"{atlas}_base_cell"
        if atlas_group_name not in root:
            raise ValueError(f"Atlas group '{atlas_group_name}' not found in zarr store")
        
        atlas_group = root[atlas_group_name]
        
        # Read gene names from zarr store
        if "var" not in atlas_group or "gene_names" not in atlas_group["var"]:
            raise ValueError("Gene names not found in zarr store. Please initialize the store first.")
        
        gene_set = list(atlas_group["var"]["gene_names"][:])
        gene_set_size = len(gene_set)
        
        print(f"Found {gene_set_size} genes in zarr store")
        
        # Get unique cell types from cellxgene
        with cellxgene_census.open_soma(census_version="2025-01-30") as census:
            cell_meta_data = cellxgene_census.get_obs(
                census, organism, column_names=["cell_type"]
            )
            if keep_all:
                cell_types = list(cell_meta_data["cell_type"].unique())
            else:
                cell_types = sorted(list(set(cell_meta_data["cell_type"])))
        
        n_cell_types = len(cell_types)
        print(f"Computing mean expression for {n_cell_types} cell types")
        
        # Check if base_cell array exists and has correct shape
        if "base_cell" not in atlas_group:
            raise ValueError("base_cell array not found in zarr store. Please initialize the store first.")
        
        base_cell_array = atlas_group["base_cell"]
        expected_shape = (n_cell_types + 1, gene_set_size)  # +1 for overall mean at row 0
        
        # Check array dimensions
        if base_cell_array.shape[1] != gene_set_size:
            raise ValueError(
                f"Gene dimension mismatch: store has {base_cell_array.shape[1]} genes, "
                f"but gene_names has {gene_set_size} genes"
            )
        
        if base_cell_array.shape[0] < expected_shape[0]:
            raise ValueError(
                f"Array too small: base_cell array has {base_cell_array.shape[0]} rows, "
                f"but need {expected_shape[0]} rows (1 for overall mean + {n_cell_types} for cell types). "
                f"Please reinitialize the store with the correct number of cell types."
            )
        
        # Compute mean expression for each cell type
        for idx, cell_type in enumerate(cell_types):
            print(f"Processing cell type {idx + 1}/{n_cell_types}: {cell_type}")
            cell_type_mean = base_cellxgene(gene_set, organism, [cell_type])
            
            # Store in zarr array (row idx + 1, since row 0 is overall mean)
            row_idx = idx + 1
            base_cell_array[row_idx, :] = cell_type_mean
        
        # Update cell_id array with cell type names
        # Note: zarr arrays can't be resized, so we delete and recreate if size differs
        if "obs" in atlas_group and "cell_id" in atlas_group["obs"]:
            cell_id_array = atlas_group["obs"]["cell_id"]
            if len(cell_id_array) != n_cell_types:
                # Delete and recreate with correct size
                del atlas_group["obs"]["cell_id"]
                gene_chunk_size, cell_chunk_size = compute_chunk_size(
                    n_cell_types, gene_set_size, max_chunk_size, min_chunk_size
                )
                atlas_group["obs"].create_array(
                    "cell_id",
                    data=cell_types,
                    shape=None,
                    dtype=None,
                    chunks=(cell_chunk_size,),
                )
            else:
                # Update in place
                cell_id_array[:] = cell_types
        else:
            # Create cell_id array
            if "obs" not in atlas_group:
                atlas_group.create_group("obs")
            gene_chunk_size, cell_chunk_size = compute_chunk_size(
                n_cell_types, gene_set_size, max_chunk_size, min_chunk_size
            )
            atlas_group["obs"].create_array(
                "cell_id",
                data=cell_types,
                shape=None,
                dtype=None,
                chunks=(cell_chunk_size,),
            )
        
        # Update metadata
        atlas_group.attrs["n_cells"] = n_cell_types
        root.attrs["n_cells"] = n_cell_types
        
        print(f"Successfully computed and stored mean expression for {n_cell_types} cell types")
        
    else:
        raise ValueError(f"Unknown atlas: {atlas}")


def _compute_euclidean_distance(x: np.ndarray, y: np.ndarray) -> float:
    """Compute Euclidean distance between two vectors."""
    return np.linalg.norm(x - y)


def _compute_cosine_distance(x: np.ndarray, y: np.ndarray) -> float:
    """Compute cosine distance between two vectors."""
    dot_product = np.dot(x, y)
    norm_x = np.linalg.norm(x)
    norm_y = np.linalg.norm(y)
    if norm_x == 0 or norm_y == 0:
        return 1.0  # Maximum distance if one vector is zero
    cosine_sim = dot_product / (norm_x * norm_y)
    return 1 - cosine_sim  # Convert similarity to distance


def _compute_manhattan_distance(x: np.ndarray, y: np.ndarray) -> float:
    """Compute Manhattan (L1) distance between two vectors."""
    return np.sum(np.abs(x - y))


def _compute_wasserstein_distance(x: np.ndarray, y: np.ndarray) -> float:
    """Compute Wasserstein distance between two vectors."""
    if wasserstein_distance is None:
        raise ImportError("scipy is required for Wasserstein distance computation")
    # For 1D distributions, we can use scipy's wasserstein_distance
    # We need to normalize to create proper probability distributions
    x_norm = x / (np.sum(x) + 1e-10)  # Add small epsilon to avoid division by zero
    y_norm = y / (np.sum(y) + 1e-10)
    # Create positions for the distributions (gene indices)
    positions = np.arange(len(x))
    return wasserstein_distance(positions, positions, x_norm, y_norm)


def cell_diff(
    zarr_store,
    atlas: str = "cellxgene",
    distance_metrics: list[Literal["euclidean", "cosine", "manhattan", "wasserstein"]] | None = None,
    max_chunk_size: int = 5000,
    min_chunk_size: int = 1000,
):
    """
    Compute the difference and distances between the overall base_cell and each cell type's base_cell.
    
    This function:
    1. Computes overall distances using multiple metrics (Euclidean, Cosine, Manhattan, Wasserstein)
    2. Stores distances in a separate 'cell_distance' group
    3. Computes normalized difference (cell_type_base - overall_base_cell) and replaces mean expression
    4. Stores normalized differences in the base_cell array (replacing original mean expression)

    Parameters:
    -----------
    zarr_store : str | zarr.Group
        Path to zarr store or zarr Group object
    atlas : str
        Atlas name (default: "cellxgene")
    distance_metrics : list[str] | None
        List of distance metrics to compute. Options: "euclidean", "cosine", "manhattan", "wasserstein".
        If None, computes all available metrics.
    max_chunk_size : int
        Maximum chunk size for zarr arrays (default: 5000)
    min_chunk_size : int
        Minimum chunk size for zarr arrays (default: 1000)

    Returns:
    --------
    None
        Results are stored in the zarr store.
    """
    # Open zarr store if path is provided
    if isinstance(zarr_store, str):
        store = zarr.storage.LocalStore(zarr_store, read_only=False)
        root = zarr.group(store)
    else:
        root = zarr_store

    if atlas == "cellxgene":
        # Get the atlas group
        atlas_group_name = f"{atlas}_base_cell"
        if atlas_group_name not in root:
            raise ValueError(f"Atlas group '{atlas_group_name}' not found in zarr store")
        
        atlas_group = root[atlas_group_name]
        
        # Check if base_cell array exists
        if "base_cell" not in atlas_group:
            raise ValueError("base_cell array not found in zarr store. Please run cell_type_base first.")
        
        base_cell_array = atlas_group["base_cell"]
        n_rows, n_genes = base_cell_array.shape
        
        if n_rows < 2:
            raise ValueError(
                "base_cell array must have at least 2 rows (1 for overall mean + at least 1 cell type). "
                "Please run cell_type_base first."
            )
        
        # Get the overall base_cell (row 0) - keep original for reference
        overall_base_cell = base_cell_array[0, :].copy()
        n_cell_types = n_rows - 1  # Subtract 1 for the overall mean row
        
        # Save original mean expressions before computing distances
        # (in case this function is called multiple times)
        original_cell_type_means = np.zeros((n_cell_types, n_genes), dtype=np.float64)
        for idx in range(n_cell_types):
            row_idx = idx + 1
            original_cell_type_means[idx, :] = base_cell_array[row_idx, :].copy()
        
        # Get cell type names
        if "obs" in atlas_group and "cell_id" in atlas_group["obs"]:
            cell_types = list(atlas_group["obs"]["cell_id"][:])
        else:
            cell_types = [f"cell_type_{i}" for i in range(n_cell_types)]
        
        print(f"Computing cell differences and distances for {n_cell_types} cell types")
        print(f"Overall base_cell shape: {overall_base_cell.shape}")
        
        # Set default distance metrics
        if distance_metrics is None:
            distance_metrics = ["euclidean", "cosine", "manhattan"]
            if wasserstein_distance is not None:
                distance_metrics.append("wasserstein")
        
        # Compute chunk sizes
        gene_chunk_size, cell_chunk_size = compute_chunk_size(
            n_cell_types, n_genes, max_chunk_size, min_chunk_size
        )
        
        # Create cell_distance group
        if "cell_distance" in root:
            distance_group = root["cell_distance"]
        else:
            distance_group = root.create_group("cell_distance")
        
        # Distance computation functions
        distance_functions = {
            "euclidean": _compute_euclidean_distance,
            "cosine": _compute_cosine_distance,
            "manhattan": _compute_manhattan_distance,
        }
        if wasserstein_distance is not None:
            distance_functions["wasserstein"] = _compute_wasserstein_distance
        
        # Compute distances for each metric using original mean expressions
        distances = {}
        for metric in distance_metrics:
            if metric not in distance_functions:
                print(f"Warning: Metric '{metric}' not available, skipping.")
                continue
            
            print(f"Computing {metric} distances...")
            dist_func = distance_functions[metric]
            metric_distances = np.zeros(n_cell_types, dtype=np.float64)
            
            for idx in range(n_cell_types):
                cell_type_base = original_cell_type_means[idx, :]
                metric_distances[idx] = dist_func(overall_base_cell, cell_type_base)
            
            # Store distances in zarr
            if metric in distance_group:
                dist_array = distance_group[metric]
                if dist_array.shape != (n_cell_types,):
                    del distance_group[metric]
                    dist_array = distance_group.create_array(
                        metric,
                        data=metric_distances,
                        shape=(n_cell_types,),
                        dtype=np.float64,
                        chunks=(cell_chunk_size,),
                    )
                else:
                    dist_array[:] = metric_distances
            else:
                distance_group.create_array(
                    metric,
                    data=metric_distances,
                    shape=(n_cell_types,),
                    dtype=np.float64,
                    chunks=(cell_chunk_size,),
                )
            
            distances[metric] = metric_distances
            print(f"  {metric} distances computed: min={metric_distances.min():.4f}, max={metric_distances.max():.4f}, mean={metric_distances.mean():.4f}")
        
        # Store cell type names in distance group
        if "cell_types" in distance_group:
            distance_group["cell_types"][:] = cell_types
        else:
            distance_group.create_array(
                "cell_types",
                data=cell_types,
                shape=None,
                dtype=None,
                chunks=(cell_chunk_size,),
            )
        
        # Compute normalized differences using original mean expressions
        # This shows if expression goes up (positive) or down (negative) compared to base
        print("Computing normalized differences...")
        
        # Normalize by the overall base_cell to get relative changes
        # Avoid division by zero
        overall_base_cell_normalized = overall_base_cell.copy()
        overall_base_cell_normalized[overall_base_cell_normalized == 0] = 1.0  # Avoid division by zero
        
        normalized_diffs = np.zeros((n_cell_types, n_genes), dtype=np.float64)
        
        for idx in range(n_cell_types):
            cell_type_base = original_cell_type_means[idx, :]
            
            # Compute normalized difference: (cell_type - base) / base
            # This gives percentage change: positive = up, negative = down
            normalized_diff = (cell_type_base - overall_base_cell) / overall_base_cell_normalized
            normalized_diffs[idx, :] = normalized_diff
        
        # Replace mean expression in base_cell array with normalized differences
        # Row 0 remains the overall base_cell, rows 1+ become normalized differences
        for idx in range(n_cell_types):
            row_idx = idx + 1
            base_cell_array[row_idx, :] = normalized_diffs[idx, :]
        
        # Store metadata
        atlas_group.attrs["cell_diff_computed"] = True
        atlas_group.attrs["normalized_differences"] = True
        distance_group.attrs["metrics"] = distance_metrics
        distance_group.attrs["n_cell_types"] = n_cell_types
        
        print(f"Successfully computed and stored:")
        print(f"  - Distances for {len(distances)} metrics")
        print(f"  - Normalized differences for {n_cell_types} cell types")
        print(f"  - Normalized differences stored in base_cell array (rows 1+)")
        
    else:
        raise ValueError(f"Unknown atlas: {atlas}")


def initialize_base_cell_store(
    gene_set,
    cell_set,
    gene_mean,
    atlas: str = "cellxgene",
    zarr_path: str = ".",
    max_chunk_size: int = 5000,
    min_chunk_size: int = 1000,
):
    if ".zarr" not in zarr_path:
        directory = zarr_path if zarr_path != "." else "."
        os.makedirs(directory, exist_ok=True)
        store_name = f"base_cell_{atlas}.zarr"
        zarr_path = os.path.join(directory, store_name)
        print('No zarr store name provided - Generic name generated')
        print(f'zarr store location: {zarr_path}')
        
    store = zarr.storage.LocalStore(zarr_path, read_only=False)
    root = zarr.create_group(store)

    cell_set_size = len(cell_set)
    gene_set_size = len(gene_set)
    cell_chunk_size, gene_chunk_size = compute_chunk_size(
        cell_set_size, gene_set_size, max_chunk_size, min_chunk_size
    )
    atlas_group = root.create_group(f'{atlas}_base_cell')
    atlas_group.create_array(
        "base_cell",
        shape=(cell_set_size + 1, gene_set_size),
        dtype=np.float64,
        chunks=(cell_chunk_size, gene_chunk_size),
        fill_value=0.0,
    )
    atlas_group["base_cell"][0,:] = gene_mean
    # either cells types or all cells but better to use cell types
    cell_id = atlas_group.create_group("obs")
    cell_id.create_array(
        "cell_id",
        data=cell_set,
        shape=None,
        dtype=None,
        chunks=(cell_chunk_size,),
    )
    # Gene names stored seperately
    gene_id = atlas_group.create_group("var")
    gene_id.create_array(
        "gene_names",
        data=gene_set,
        shape=None,
        dtype=None,
        chunks=(gene_chunk_size,),
    )
    atlas_group.attrs["n_cells"] = cell_set_size
    atlas_group.attrs["n_genes"] = gene_set_size
    root.attrs["n_cells"] = cell_set_size
    root.attrs["n_genes"] = gene_set_size

    return root, zarr_path


def plot_heatmap(
    zarr_store,
    atlas: str = "cellxgene",
    output_path: str | None = None,
    max_genes: int | None = None,
    max_cell_types: int | None = None,
):
    """
    Plot a heatmap of normalized differences from the zarr store.
    
    The heatmap shows normalized differences (up/down regulation) for each cell type
    compared to the base cell. Positive values (red) indicate up-regulation,
    negative values (blue) indicate down-regulation.

    Parameters:
    -----------
    zarr_store : str | zarr.Group
        Path to zarr store or zarr Group object
    atlas : str
        Atlas name (default: "cellxgene")
    output_path : str | None
        Path to save the plot. If None, displays the plot.
    max_genes : int | None
        Maximum number of genes to display. If None, displays all.
    max_cell_types : int | None
        Maximum number of cell types to display. If None, displays all.

    Returns:
    --------
    None
        Displays or saves the heatmap plot.
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        raise ImportError("matplotlib and seaborn are required for plotting. Install with: pip install matplotlib seaborn")
    
    # Open zarr store if path is provided
    if isinstance(zarr_store, str):
        store = zarr.storage.LocalStore(zarr_store, read_only=False)
        root = zarr.group(store)
    else:
        root = zarr_store

    if atlas == "cellxgene":
        # Get the atlas group
        atlas_group_name = f"{atlas}_base_cell"
        if atlas_group_name not in root:
            raise ValueError(f"Atlas group '{atlas_group_name}' not found in zarr store")
        
        atlas_group = root[atlas_group_name]
        
        # Check if base_cell array exists
        if "base_cell" not in atlas_group:
            raise ValueError("base_cell array not found in zarr store.")
        
        base_cell_array = atlas_group["base_cell"]
        n_rows, n_genes = base_cell_array.shape
        
        if n_rows < 2:
            raise ValueError("base_cell array must have at least 2 rows.")
        
        # Get normalized differences (rows 1+)
        n_cell_types = n_rows - 1
        normalized_diffs = base_cell_array[1:, :]  # Skip row 0 (overall base)
        
        # Get gene names
        if "var" in atlas_group and "gene_names" in atlas_group["var"]:
            gene_names = list(atlas_group["var"]["gene_names"][:])
        else:
            gene_names = [f"gene_{i}" for i in range(n_genes)]
        
        # Get cell type names
        if "obs" in atlas_group and "cell_id" in atlas_group["obs"]:
            cell_types = list(atlas_group["obs"]["cell_id"][:])
        else:
            cell_types = [f"cell_type_{i}" for i in range(n_cell_types)]
        
        # Limit display size if requested
        if max_genes is not None and n_genes > max_genes:
            # Select top varying genes
            gene_variance = np.var(normalized_diffs, axis=0)
            top_gene_indices = np.argsort(gene_variance)[-max_genes:]
            normalized_diffs = normalized_diffs[:, top_gene_indices]
            gene_names = [gene_names[i] for i in top_gene_indices]
            print(f"Displaying top {max_genes} most varying genes")
        
        if max_cell_types is not None and n_cell_types > max_cell_types:
            normalized_diffs = normalized_diffs[:max_cell_types, :]
            cell_types = cell_types[:max_cell_types]
            print(f"Displaying first {max_cell_types} cell types")
        
        # Create heatmap
        plt.figure(figsize=(max(12, len(gene_names) * 0.1), max(8, len(cell_types) * 0.3)))
        
        # Use a diverging colormap: blue (negative/down) -> white (zero) -> red (positive/up)
        vmax = np.max(np.abs(normalized_diffs))
        vmin = -vmax
        
        sns.heatmap(
            normalized_diffs,
            xticklabels=gene_names,
            yticklabels=cell_types,
            cmap="RdBu_r",  # Red-Blue reversed: red=up, blue=down
            center=0,
            vmin=vmin,
            vmax=vmax,
            cbar_kws={"label": "Normalized Difference (Up/Down Regulation)"},
            linewidths=0.5,
            linecolor="gray",
        )
        
        plt.title("Cell Type Expression Differences (Normalized)\nRed = Up-regulation, Blue = Down-regulation", 
                  fontsize=14, pad=20)
        plt.xlabel("Genes", fontsize=12)
        plt.ylabel("Cell Types", fontsize=12)
        plt.xticks(rotation=90, ha="right", fontsize=8)
        plt.yticks(rotation=0, fontsize=8)
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches="tight")
            print(f"Heatmap saved to {output_path}")
        else:
            plt.show()
        
        plt.close()
        
    else:
        raise ValueError(f"Unknown atlas: {atlas}")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Cell type base expression and difference computation"
    )
    parser.add_argument(
        "command",
        choices=["init", "cell_type_base", "cell_diff", "plot"],
        help="Command to execute"
    )
    parser.add_argument(
        "--zarr-store",
        type=str,
        required=True,
        help="Path to zarr store"
    )
    parser.add_argument(
        "--atlas",
        type=str,
        default="cellxgene",
        help="Atlas name (default: cellxgene)"
    )
    parser.add_argument(
        "--organism",
        type=str,
        default="mus_musculus",
        help="Organism name (default: mus_musculus)"
    )
    parser.add_argument(
        "--keep-all",
        action="store_true",
        help="Keep all cell types (for cell_type_base)"
    )
    parser.add_argument(
        "--max-chunk-size",
        type=int,
        default=5000,
        help="Maximum chunk size for zarr arrays (default: 5000)"
    )
    parser.add_argument(
        "--min-chunk-size",
        type=int,
        default=1000,
        help="Minimum chunk size for zarr arrays (default: 1000)"
    )
    parser.add_argument(
        "--distance-metrics",
        nargs="+",
        choices=["euclidean", "cosine", "manhattan", "wasserstein"],
        help="Distance metrics to compute (for cell_diff). If not specified, computes all available."
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output path for plot (for plot command)"
    )
    parser.add_argument(
        "--max-genes",
        type=int,
        help="Maximum number of genes to display in heatmap (for plot command)"
    )
    parser.add_argument(
        "--max-cell-types",
        type=int,
        help="Maximum number of cell types to display in heatmap (for plot command)"
    )
    
    args = parser.parse_args()
    
    if args.command == "init":
        print("Initializing base cell store...")
        gene_set, base_cell_mean = base_cell(
            atlas=args.atlas,
            organism=args.organism,
        )
        
        # Get cell types for initialization
        import cellxgene_census
        with cellxgene_census.open_soma(census_version="2025-01-30") as census:
            cell_meta_data = cellxgene_census.get_obs(
                census, args.organism, column_names=["cell_type"]
            )
            cell_types = sorted(list(set(cell_meta_data["cell_type"])))
        
        root, zarr_path = initialize_base_cell_store(
            gene_set=list(gene_set),
            cell_set=cell_types,
            gene_mean=base_cell_mean,
            atlas=args.atlas,
            zarr_path=args.zarr_store,
            max_chunk_size=args.max_chunk_size,
            min_chunk_size=args.min_chunk_size,
        )
        print(f"Base cell store initialized at: {zarr_path}")
    
    elif args.command == "cell_type_base":
        print("Computing cell type base expressions...")
        cell_type_base(
            zarr_store=args.zarr_store,
            atlas=args.atlas,
            organism=args.organism,
            keep_all=args.keep_all,
            max_chunk_size=args.max_chunk_size,
            min_chunk_size=args.min_chunk_size,
        )
        print("Cell type base expressions computed.")
    
    elif args.command == "cell_diff":
        print("Computing cell differences and distances...")
        cell_diff(
            zarr_store=args.zarr_store,
            atlas=args.atlas,
            distance_metrics=args.distance_metrics,
            max_chunk_size=args.max_chunk_size,
            min_chunk_size=args.min_chunk_size,
        )
        print("Cell differences and distances computed.")
    
    elif args.command == "plot":
        print("Generating heatmap...")
        plot_heatmap(
            zarr_store=args.zarr_store,
            atlas=args.atlas,
            output_path=args.output,
            max_genes=args.max_genes,
            max_cell_types=args.max_cell_types,
        )
        print("Heatmap generated.")
