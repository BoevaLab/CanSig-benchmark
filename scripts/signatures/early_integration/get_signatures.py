import logging
import pathlib as pl
from argparse import ArgumentParser

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from utils import read_adata, get_args, set_seed


logging.basicConfig()
logging.root.setLevel(logging.INFO)

_LOGGER = logging.getLogger(__name__)
_OFFSET: int = 10_000
_LATENT = "latent"

_ARGS: list = [
    ("-n", "--n-cluster", int),
    ("-l", "--latent", str)
]


def add_latent(path: str, adata: sc.AnnData) -> sc.AnnData:
    if path.endswith(".npz"):
        _LOGGER.info(f"Loading graph from {path}")
        data = np.load(path, allow_pickle=True)
        # Numpy saves everything as arrays. .item retrieves the original object.
        adata = adata[data["index"]].copy()
        adata.obsp[f"connectivities"] = data["connectivities"].item()
        adata.obsp[f"distances"] = data["distances"].item()
        adata.uns["neighbors"] = {'connectivities_key': f'connectivities',
                                       'distances_key': f'distances',
                                       'params': {'n_neighbors': 15, 'method': 'umap', 'random_state': 0,
                                                  'metric': 'euclidean'}}
    elif path.endswith(".csv"):
        _LOGGER.info(f"Loading latents from {path}")

        latent = pd.read_csv(path, index_col=0)
        adata = adata[latent.index].copy()
        if np.any(latent.index != adata.obs_names):
            raise ValueError("Mismatch in index of the latent space and the adata. Something went wrong.")
        adata.obsm[_LATENT] = latent.values
        _LOGGER.info(f"Computing knn-graphs")
        sc.pp.neighbors(adata, use_rep=_LATENT, n_neighbors=15)
    else:
        raise ValueError(f"Unknown file type, {path}.")
    return adata



def get_cluster(adata, n_cluster, random_seed ,start: float = 1e-4, end: float= 2., epsilon: float = 1e-8):
    for i in range(10):
        try:
            while end - start > epsilon:
                mid = (end + start) / 2.0
                sc.tl.leiden(adata, resolution=mid, random_state=_OFFSET*i+random_seed)
                n_tmp = adata.obs["leiden"].nunique()
                if n_tmp == n_cluster:
                    _LOGGER.info(f"Found {n_tmp} cluster.")
                    break
                if n_tmp > n_cluster:
                    end = mid
                if n_tmp < n_cluster:
                    start = mid
            else:
                raise ValueError("Number of clusters doesn't match.")
        except ValueError:
            _LOGGER.info(f"Unsuccessful for random state {i}. Trying next random state.")
            pass
        else:
            break



def get_metasigs(adata: ad.AnnData, latent: str, n_cluster: int, random_seed: int) -> pd.DataFrame:

    CLUSTERS = "leiden"
    NAMES = "names"
    QVALUE = "pvals_adj"
    LOGFOLD = "logfoldchanges"
    _LATENT = "latent"
    _LOGGER.info(f"Running clustering with {n_cluster} clusters.")

    adata = add_latent(latent, adata)

    get_cluster(adata, n_cluster=n_cluster, random_seed=random_seed)

    _LOGGER.info("Running DE analysis.")
    n_cell_per_cluster = adata.obs[CLUSTERS].value_counts()
    cluster_for_de = n_cell_per_cluster[n_cell_per_cluster>5].index.to_list()
    sc.tl.rank_genes_groups(adata, groupby=CLUSTERS, groups=cluster_for_de, n_genes=500)

    _LOGGER.info("Getting meta-signatures.")
    metasigs = []
    for n_metasig, cluster in enumerate(cluster_for_de):
        _LOGGER.info(f" Running DE for {cluster}.")
        gene_list = sc.get.rank_genes_groups_df(adata, group=str(cluster))
        gene_list = gene_list[(gene_list[LOGFOLD] > 0) & (gene_list[QVALUE] < 0.05)]
        gene_list = gene_list.head(n=50)[NAMES]
        gene_list.name = f"Sig. {n_metasig+1}"
        metasigs.append(gene_list[:50])

    return pd.concat(metasigs, axis=1)


def main():

    args = get_args(_ARGS)
    set_seed(args.random_seed)
    adata = read_adata(args.input, n_hvg=None)
    _LOGGER.info(f"Loading latent codes from {args.latent}")
    metasigs = get_metasigs(adata, args.latent, n_cluster=args.n_cluster, random_seed=args.random_seed)
    metasigs.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
