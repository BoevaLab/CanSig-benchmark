import argparse
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

logger = logging.getLogger(__name__)

_PHASE = "phase"
_G1 = "G1"
_LATENT = "_cansig_latent"
_LEIDEN = "_cansig_leiden"
_CELLTYPE = "celltype"
_SAMPLE = "sample"
_15_NEIGHBORS_KEY = "15"
_50_NEIGHBORS_KEY = "50"
_90_NEIGHBORS_KEY = "90"



def convert_knn_graph_to_idx(X: csr_matrix) -> tuple[np.ndarray, np.ndarray]:
    """Convert a kNN graph to indices and distances."""
    check_array(X, accept_sparse="csr")
    check_square(X)

    n_neighbors = np.unique(X.nonzero()[0], return_counts=True)[1]
    n_neighbors = int(np.min((np.unique(n_neighbors))))
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Precomputed sparse input")
        nn_obj = NearestNeighbors(n_neighbors=n_neighbors, metric="precomputed").fit(X)
        kneighbors = nn_obj.kneighbors(X)
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
    labels = np.asarray(pd.Categorical(labels).codes)
    knn_dists, knn_idx = convert_knn_graph_to_idx(X)

    if perplexity is None:
        perplexity = np.floor(knn_idx.shape[1] / 3)

    n_labels = len(np.unique(labels))

    simpson = compute_simpson_index(knn_dists, knn_idx, labels, n_labels, perplexity=perplexity)
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
    batches = np.asarray(pd.Categorical(batches).codes)
    lisi = lisi_knn(X, batches, perplexity=perplexity)
    ilisi = np.nanmedian(lisi)
    if scale:
        nbatches = len(np.unique(batches))
        ilisi = (ilisi - 1) / (nbatches - 1)
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
    labels = np.asarray(pd.Categorical(labels).codes)
    lisi = lisi_knn(X, labels, perplexity=perplexity)
    clisi = np.nanmedian(lisi)
    if scale:
        nlabels = len(np.unique(labels))
        clisi = (nlabels - clisi) / (nlabels - 1)
    return clisi


def set_neighbors(adata: sc.AnnData, distances, connectivities: np.ndarray, n_neighbors: int) -> sc.AnnData:
    logger.info(f"Setting {n_neighbors}-nn Graph.")
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

    return distances, connectivites


def add_latent(path: str, adata: sc.AnnData) -> sc.AnnData:
    if path.endswith(".npz"):
        logger.info(f"Loading graph from {path}")
        data = np.load(path, allow_pickle=True)
        # Numpy saves everything as arrays. .item retrieves the original object.
        adata = set_neighbors(adata, connectivities=data["connectivities"].item(), distances=data["distances"].item(), n_neighbors=15)
        distances, connectivities = diffusion_nn(adata.obsp["15_connectivities"], 50)
        adata = set_neighbors(adata, connectivities=connectivities, distances=distances, n_neighbors=50)
        distances, connectivities = diffusion_nn(adata.obsp["15_connectivities"], 90)
        adata = set_neighbors(adata, connectivities=connectivities, distances=distances, n_neighbors=90)
    elif path.endswith(".csv"):
        logger.info(f"Loading latents from {path}")

        latent = pd.read_csv(path, index_col=0)

        if len(latent.index) != len(adata.obs_names) or set(adata.obs_names) != set(latent.index):
            raise ValueError("Mismatch in index of the latent space and the adata. Something went wrong.")

        latent = latent.loc[adata.obs_names]
        if np.any(latent.index != adata.obs_names):
            raise ValueError("Mismatch in index of the latent space and the adata. Something went wrong.")
        adata.obsm[_LATENT] = latent.values
        logger.info(f"Computing knn-graphs")
        sc.pp.neighbors(adata, use_rep=_LATENT, n_neighbors=15, key_added="15")
        sc.pp.neighbors(adata, use_rep=_LATENT, n_neighbors=50, key_added="50")
        sc.pp.neighbors(adata, use_rep=_LATENT, n_neighbors=90, key_added="90")
    else:
        raise ValueError(f"Unknown file type, {path}.")
    return adata


def get_ari_nmi(adata: sc.AnnData) -> pd.DataFrame:
    logger.info("Computing ari and nmi")
    metrics = {"ARI": adjusted_rand_score,
               "NMI": normalized_mutual_info_score}
    resolutions = [0.05 * i for i in range(2, 21)]
    n_randoms = 5
    results = pd.DataFrame(None, columns=list(range(n_randoms)),
                           index=pd.MultiIndex.from_product([metrics.keys(), resolutions]))

    for random_seed in range(n_randoms):
        for n_res, res in enumerate(resolutions):
            sc.tl.leiden(adata, resolution=res, random_state=random_seed, key_added=_LEIDEN,
                         neighbors_key=_15_NEIGHBORS_KEY)
            for metric_name, metric in metrics.items():
                results.loc[(metric_name, res), random_seed] = metric(adata.obs[_LEIDEN],
                                                                      adata.obs[_CELLTYPE])

    return results


def get_kbets(adata) -> dict:

    labels = adata.obs[_CELLTYPE].values
    batchs = adata.obs[_SAMPLE].values
    results = {"kbet_per_label": kbet_per_label(adata.obsp[_50_NEIGHBORS_KEY + "_connectivities"],
                                               batchs, labels)}
    return results


def get_clisi_knn(adata: anndata.AnnData) -> dict:
    clisi = clisi_knn(adata.obsp[_90_NEIGHBORS_KEY + "_distances"], adata.obs[_CELLTYPE].values)
    return {"clisi": clisi}
def get_ilisis_knn(adata: anndata.AnnData) -> dict:
    ilisi = ilisi_knn(adata.obsp[_90_NEIGHBORS_KEY + "_distances"], adata.obs[_SAMPLE].values)
    return {"ilisis": ilisi}

def get_graph_connectivity(adata: anndata.AnnData) -> dict:
    graph_connect = graph_connectivity(adata.obsp[_15_NEIGHBORS_KEY + "_distances"], adata.obs[_CELLTYPE].values)
    return {"graph_connect": graph_connect}

def save_results(output_cluster: str, output_metrics, results: dict, clustering_results: pd.DataFrame) -> None:
    pathlib.Path(output_cluster).parent.mkdir(exist_ok=True, parents=True)
    clustering_results.to_csv(output_cluster)
    pd.DataFrame(results, index=[0]).to_csv(output_metrics)


def compute_umap(adata) -> pd.DataFrame:
    sc.tl.umap(adata, neighbors_key="15")
    df = pd.DataFrame(adata.obsm["X_umap"], columns=[f"UMAP_{i}" for i in range(2)])
    df.index = adata.obs_names
    return df



def save_umaps(df: pd.DataFrame, path: str) -> None:
    df.to_csv(path)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-path", type=str)
    parser.add_argument("--latent", type=str)
    parser.add_argument("--output-cluster", type=str)
    parser.add_argument("--output-metrics", type=str)
    parser.add_argument("--output-umap", type=str)
    parser.add_argument("--plotting", type=bool)

    return parser.parse_args()


def main():
    args = get_args()
    adata = read_adata(args.data_path, n_hvg=None, batch_key=False, remove_cycling=True)
    adata = add_latent(args.latent, adata)
    results = get_kbets(adata)
    results.update(get_ilisis_knn(adata))
    results.update(get_clisi_knn(adata))
    results.update(get_graph_connectivity(adata))
    clustering_results = get_ari_nmi(adata)

    save_results(args.output_cluster, args.output_metrics, results, clustering_results)

    if args.plotting:
        umap = compute_umap(adata)
        save_umaps(umap, args.output_umap)
    else:
        pd.DataFrame().to_csv(args.output_umap)




if __name__ == '__main__':
    main()
