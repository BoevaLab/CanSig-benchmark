import argparse
from itertools import product
import logging
import pathlib
import warnings

import anndata
import numpy as np
import pandas as pd
import pynndescent
import scanpy as sc
from scanpy.neighbors import _compute_connectivities_umap
from scib_metrics import graph_connectivity, clisi_knn, ilisi_knn
from scib_metrics import kbet_per_label
from scib_metrics.utils import check_square, compute_simpson_index
from scib_metrics.utils._diffusion_nn import _compute_transitions, _compute_eigen
from scipy.sparse import csr_matrix
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.neighbors import NearestNeighbors
from sklearn.utils import check_array

from utils import read_adata


import warnings
from numba.core.errors import NumbaDeprecationWarning
warnings.filterwarnings("ignore", category=NumbaDeprecationWarning)


logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

_LOGGER = logging.getLogger(__name__)

_LATENT = "_cansig_latent"
_LEIDEN = "_cansig_leiden"
_CELLTYPE = "celltype"
_SAMPLE = "sample"
_15_NEIGHBORS_KEY = "15"
_50_NEIGHBORS_KEY = "50"
_90_NEIGHBORS_KEY = "90"



def convert_knn_graph_to_idx(X: csr_matrix) -> tuple[np.ndarray, np.ndarray]:
    """Convert a kNN graph to indices and distances."""
    _LOGGER.info("Converting kNN graph to indices and distances")
    check_array(X, accept_sparse="csr")
    check_square(X)

    n_neighbors = np.unique(X.nonzero()[0], return_counts=True)[1]
    n_neighbors = int(np.min((np.unique(n_neighbors))))
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Precomputed sparse input")
        nn_obj = NearestNeighbors(n_neighbors=n_neighbors, metric="precomputed").fit(X)
        kneighbors = nn_obj.kneighbors(X)
    _LOGGER.info(f"Converted kNN graph with {n_neighbors} neighbors")
    return kneighbors


def lisi_knn(X: csr_matrix, labels: np.ndarray, perplexity: float = None) -> np.ndarray:
    """Compute the local inverse simpson index (LISI) for each cell :cite:p:`korsunsky2019harmony`.

    Parameters
    ----------
    X
        Array of shape (n_cells, n_cells) with non-zero values
        representing distances to exactly each cell's k nearest neighbors.
    labels
        Array of shape (n_cells,) representing label values
        for each cell.
    perplexity
        Parameter controlling effective neighborhood size. If None, the
        perplexity is set to the number of neighbors // 3.

    Returns
    -------
    lisi
        Array of shape (n_cells,) with the LISI score for each cell.
    """
    _LOGGER.info("Computing local inverse simpson index (LISI)")
    labels = np.asarray(pd.Categorical(labels).codes)
    knn_dists, knn_idx = convert_knn_graph_to_idx(X)

    if perplexity is None:
        perplexity = np.floor(knn_idx.shape[1] / 3)

    n_labels = len(np.unique(labels))

    simpson = compute_simpson_index(knn_dists, knn_idx, labels, n_labels, perplexity=perplexity)
    _LOGGER.info(f"Computed LISI for {len(labels)} cells with {n_labels} unique labels")
    return 1 / simpson


def ilisi_knn(X: csr_matrix, batches: np.ndarray, perplexity: float = None, scale: bool = True) -> np.ndarray:
    """Compute the integration local inverse simpson index (iLISI) for each cell :cite:p:`korsunsky2019harmony`.

    Returns a scaled version of the iLISI score for each cell, by default :cite:p:`luecken2022benchmarking`.

    Parameters
    ----------
    X
        Array of shape (n_cells, n_cells) with non-zero values
        representing distances to exactly each cell's k nearest neighbors.
    batches
        Array of shape (n_cells,) representing batch values
        for each cell.
    perplexity
        Parameter controlling effective neighborhood size. If None, the
        perplexity is set to the number of neighbors // 3.
    scale
        Scale lisi into the range [0, 1]. If True, higher values are better.

    Returns
    -------
    ilisi
        Array of shape (n_cells,) with the iLISI score for each cell.
    """
    _LOGGER.info("Computing integration local inverse simpson index (iLISI)")
    batches = np.asarray(pd.Categorical(batches).codes)
    lisi = lisi_knn(X, batches, perplexity=perplexity)
    ilisi = np.nanmedian(lisi)
    if scale:
        nbatches = len(np.unique(batches))
        ilisi = (ilisi - 1) / (nbatches - 1)
    _LOGGER.info(f"Computed iLISI score: {ilisi:.4f}")
    return ilisi


def clisi_knn(X: csr_matrix, labels: np.ndarray, perplexity: float = None, scale: bool = True) -> np.ndarray:
    """Compute the cell-type local inverse simpson index (cLISI) for each cell :cite:p:`korsunsky2019harmony`.

    Returns a scaled version of the cLISI score for each cell, by default :cite:p:`luecken2022benchmarking`.

    Parameters
    ----------
    X
        Array of shape (n_cells, n_cells) with non-zero values
        representing distances to exactly each cell's k nearest neighbors.
    labels
        Array of shape (n_cells,) representing cell type label values
        for each cell.
    perplexity
        Parameter controlling effective neighborhood size. If None, the
        perplexity is set to the number of neighbors // 3.
    scale
        Scale lisi into the range [0, 1]. If True, higher values are better.

    Returns
    -------
    clisi
        Array of shape (n_cells,) with the cLISI score for each cell.
    """
    _LOGGER.info("Computing cell-type local inverse simpson index (cLISI)")
    labels = np.asarray(pd.Categorical(labels).codes)
    lisi = lisi_knn(X, labels, perplexity=perplexity)
    clisi = np.nanmedian(lisi)
    if scale:
        nlabels = len(np.unique(labels))
        clisi = (nlabels - clisi) / (nlabels - 1)
    _LOGGER.info(f"Computed cLISI score: {clisi:.4f}")
    return clisi


def set_neighbors(adata: sc.AnnData, distances, connectivities: np.ndarray, n_neighbors: int) -> sc.AnnData:
    _LOGGER.info(f"Setting {n_neighbors}-nn Graph.")
    adata.obsp[f"{n_neighbors}_connectivities"] = connectivities
    adata.obsp[f"{n_neighbors}_distances"] = distances

    adata.uns[f"{n_neighbors}"] = {'connectivities_key': f'{n_neighbors}_connectivities',
                                   'distances_key': f'{n_neighbors}_distances',
                                   'params': {'n_neighbors': n_neighbors, 'method': 'umap', 'random_state': 0,
                                              'metric': 'euclidean'}}

    return adata


def diffusion_nn(X: csr_matrix, k: int, n_comps: int = 100):
    """Diffusion-based neighbors.

    This function generates a nearest neighbour list from a connectivities matrix.
    This allows us to select a consistent number of nearest neighbors across all methods.

    This differs from the original scIB implemenation by leveraging diffusion maps. Here we
    embed the data with diffusion maps in which euclidean distance represents well the diffusion
    distance. We then use pynndescent to find the nearest neighbours in this embedding space.

    Parameters
    ----------
    X
        Array of shape (n_cells, n_cells) with non-zero values
        representing connectivities.
    k
        Number of nearest neighbours to select.
    n_comps
        Number of components for diffusion map

    Returns
    -------
    Neighbors graph
    """
    _LOGGER.info(f"Computing diffusion-based neighbors with k={k} and n_comps={n_comps}")
    transitions = _compute_transitions(X)
    evals, evecs = _compute_eigen(transitions, n_comps=n_comps)
    evals += 1e-8  # Avoid division by zero
    # Multiscale such that the number of steps t gets "integrated out"
    embedding = evecs
    scaled_evals = np.array([e if e == 1 else e / (1 - e) for e in evals])
    embedding *= scaled_evals
    nn_obj = pynndescent.NNDescent(embedding, n_neighbors=k + 1)
    neigh_inds, neigh_distances = nn_obj.neighbor_graph

    distances, connectivites = _compute_connectivities_umap(neigh_inds, neigh_distances, X.shape[0], k)
    _LOGGER.info("Completed diffusion-based neighbor computation")

    return distances, connectivites


def add_latent(path: str, adata: sc.AnnData) -> sc.AnnData:
    _LOGGER.info(f"Adding latent representation from {path}")
    if path.endswith(".npz"):
        _LOGGER.info(f"Loading graph from {path}")
        data = np.load(path, allow_pickle=True)
        # Numpy saves everything as arrays. .item retrieves the original object.
        adata = set_neighbors(adata, connectivities=data["connectivities"].item(), distances=data["distances"].item(), n_neighbors=15)
        distances, connectivities = diffusion_nn(adata.obsp["15_connectivities"], 50)
        adata = set_neighbors(adata, connectivities=connectivities, distances=distances, n_neighbors=50)
        distances, connectivities = diffusion_nn(adata.obsp["15_connectivities"], 90)
        adata = set_neighbors(adata, connectivities=connectivities, distances=distances, n_neighbors=90)
    elif path.endswith(".csv"):
        _LOGGER.info(f"Loading latents from {path}")

        latent = pd.read_csv(path, index_col=0)

        if len(latent.index) != len(adata.obs_names) or set(adata.obs_names) != set(latent.index):
            raise ValueError("Mismatch in index of the latent space and the adata. Something went wrong.")

        latent = latent.loc[adata.obs_names]
        if np.any(latent.index != adata.obs_names):
            raise ValueError("Mismatch in index of the latent space and the adata. Something went wrong.")
        adata.obsm[_LATENT] = latent.values
        _LOGGER.info(f"Computing knn-graphs")
        for n_neighbors in [15, 50, 90]:
            try: 
                sc.pp.neighbors(adata, use_rep=_LATENT, n_neighbors=n_neighbors, key_added=str(n_neighbors), random_state=42)
            except Exception as e:
                print(e)
                sc.pp.neighbors(adata, use_rep=_LATENT, n_neighbors=n_neighbors, key_added=str(n_neighbors), random_state=52)
                
    else:
        raise ValueError(f"Unknown file type, {path}.")
    _LOGGER.info("Completed adding latent representation")
    return adata


def get_ari_nmi(adata: sc.AnnData) -> pd.DataFrame:
    _LOGGER.info("Computing ari and nmi")
    metrics = {"ARI": adjusted_rand_score,
               "NMI": normalized_mutual_info_score}
    resolutions = [0.1 * i for i in range(2, 11)]
    n_randoms = 5
    results = pd.DataFrame(None, columns=list(range(n_randoms)),
                           index=pd.MultiIndex.from_product([metrics.keys(), resolutions]))

    for random_seed, res in product(range(n_randoms), resolutions):
        _LOGGER.info(f"Computing leiden for {res} resolution and random seed {random_seed}.")
        sc.tl.leiden(adata, resolution=res, random_state=random_seed, key_added=_LEIDEN,
                        neighbors_key=_15_NEIGHBORS_KEY)
        for metric_name, metric in metrics.items():
            metric_value = metric(adata.obs[_LEIDEN], adata.obs[_CELLTYPE])
            _LOGGER.info(f"{metric_name}: {metric_value}.")
            results.loc[(metric_name, res), random_seed] = metric_value

    return results


def get_kbets(adata) -> dict:
    _LOGGER.info("Computing k-bet per label")
    labels = adata.obs[_CELLTYPE].values
    batchs = adata.obs[_SAMPLE].values
    results = {"kbet_per_label": kbet_per_label(adata.obsp[_50_NEIGHBORS_KEY + "_connectivities"],
                                               batchs, labels)}
    _LOGGER.info("Completed k-bet computation")
    return results


def get_clisi_knn(adata: anndata.AnnData) -> dict:
    _LOGGER.info("Computing cLISI using k-NN")
    clisi = clisi_knn(adata.obsp[_90_NEIGHBORS_KEY + "_distances"], adata.obs[_CELLTYPE].values)
    _LOGGER.info("Completed cLISI computation")
    return {"clisi": clisi}


def get_ilisis_knn(adata: anndata.AnnData) -> dict:
    _LOGGER.info("Computing iLISI using k-NN")
    ilisi = ilisi_knn(adata.obsp[_90_NEIGHBORS_KEY + "_distances"], adata.obs[_SAMPLE].values)
    _LOGGER.info("Completed iLISI computation")
    return {"ilisis": ilisi}


def get_graph_connectivity(adata: anndata.AnnData) -> dict:
    _LOGGER.info("Computing graph connectivity")
    graph_connect = graph_connectivity(adata.obsp[_15_NEIGHBORS_KEY + "_distances"], adata.obs[_CELLTYPE].values)
    _LOGGER.info("Completed graph connectivity computation")
    return {"graph_connect": graph_connect}


def save_results(output, results: dict) -> None:
    _LOGGER.info(f"Saving results to {output}")
    pathlib.Path(output).parent.mkdir(exist_ok=True, parents=True)
    pd.DataFrame(results, index=[0]).to_csv(output)
    _LOGGER.info("Results saved successfully")


def compute_umap(adata) -> pd.DataFrame:
    _LOGGER.info("Computing UMAP embedding")
    sc.tl.umap(adata, neighbors_key="15")
    df = pd.DataFrame(adata.obsm["X_umap"], columns=[f"UMAP_{i}" for i in range(2)])
    df.index = adata.obs_names
    _LOGGER.info("UMAP computation completed")
    return df


def save_umaps(df: pd.DataFrame, path: str) -> None:
    _LOGGER.info(f"Saving UMAP coordinates to {path}")
    df.to_csv(path)
    _LOGGER.info("UMAP coordinates saved successfully")


def get_args():
    _LOGGER.info("Parsing command line arguments")
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-path", type=str)
    parser.add_argument("--latent", type=str)
    parser.add_argument("--output", type=str)
    parser.add_argument("--output-umap", type=str)
    parser.add_argument("--plotting", type=bool)

    return parser.parse_args()


def main():
    _LOGGER.info("Starting main execution")
    args = get_args()
    _LOGGER.info(f"Loading data from {args.data_path}")
    adata = read_adata(args.data_path, n_hvg=None, batch_key=False, remove_cycling=True)
    adata = add_latent(args.latent, adata)
    
    _LOGGER.info("Computing metrics")
    results = get_kbets(adata)
    results.update(get_ilisis_knn(adata))
    results.update(get_clisi_knn(adata))
    results.update(get_graph_connectivity(adata))
    clustering_results = get_ari_nmi(adata)
    best_clustering = clustering_results.xs(clustering_results.mean(1).groupby(level=1).mean().idxmax(), level=1)
    results.update(dict(best_clustering.mean(1)))
    save_results(args.output, results)

    if args.plotting:
        _LOGGER.info("Plotting enabled, computing UMAP")
        umap = compute_umap(adata)
        save_umaps(umap, args.output_umap)
    else:
        _LOGGER.info("Plotting disabled, saving empty DataFrame")
        pd.DataFrame().to_csv(args.output_umap)

    _LOGGER.info("Main execution completed successfully")


if __name__ == '__main__':
    main()