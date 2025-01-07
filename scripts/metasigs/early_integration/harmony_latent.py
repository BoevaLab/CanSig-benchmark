import logging
import pathlib as pl
from argparse import ArgumentParser

import anndata
import pandas as pd
import scanpy as sc


_LATENT = "latent"

logging.basicConfig()
logging.root.setLevel(logging.INFO)

_LOGGER = logging.getLogger(__name__)


def read_anndata(path: str):
    adata = sc.read(path)
    return adata


def get_latent(adata: anndata.AnnData) -> pd.DataFrame:
    BATCH_KEY = "sample"
    adata = adata[adata.obs["phase"] == "G1"].copy()
    sc.pp.highly_variable_genes(adata, n_top_genes=4000, subset=True)
    sc.tl.pca(adata)
    sc.external.pp.harmony_integrate(adata, key=BATCH_KEY, basis="X_pca", adjusted_basis=_LATENT)
    latent = pd.DataFrame(adata.obsm[_LATENT], index=adata.obs_names)
    return latent


def get_args():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True)
    parser.add_argument("-o", "--output", type=str, required=True)
    return parser.parse_args()


def main():
    args = get_args()
    output = pl.Path(args.output)
    _LOGGER.info(f"Making output dir: {str(output.parents[0])}")
    output.parents[0].mkdir(exist_ok=True, parents=True)
    _LOGGER.info(f"Reading adata from {args.input}")
    adata = read_anndata(args.input)
    latent = get_latent(adata)
    _LOGGER.info(f"Writing latent codes to {output}")
    latent.to_csv(output)


if __name__ == '__main__':
    main()
