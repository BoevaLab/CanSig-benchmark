import logging
import pathlib as pl
from argparse import ArgumentParser

import anndata
import pandas as pd
import scanpy as sc
import numpy as np

from utils import read_adata, get_args, set_seed

BATCH_KEY = "sample"

logging.basicConfig()
logging.root.setLevel(logging.INFO)

_LOGGER = logging.getLogger(__name__)


_ARGS = [("-n", "--n-hvg", str),
         ("-b", "--batch-key", str)]

def get_latent(adata: anndata.AnnData) -> pd.DataFrame | dict:
    sc.tl.pca(adata)
    sc.external.pp.bbknn(adata, batch_key=BATCH_KEY)
    latent = {"connectivities": adata.obsp["connectivities"], "distances": adata.obsp["distances"], "index": adata.obs_names}
    return latent


def main():
    args = get_args(args=_ARGS)
    set_seed(args.random_seed)
    adata = read_adata(args.input, n_hvg=args.n_hvg, batch_key=args.batch_key)
    latent = get_latent(adata)
    _LOGGER.info(f"Writing latent codes to {args.output}")
    np.savez(args.output.with_suffix(""), **latent)


if __name__ == '__main__':
    main()
