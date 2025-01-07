import logging
import pathlib as pl
from argparse import ArgumentParser

import pandas as pd
import scanpy as sc

logging.basicConfig()
logging.root.setLevel(logging.INFO)

_LOGGER = logging.getLogger(__name__)


def read_anndata(path: str):
    adata = sc.read(path)
    return adata


def get_args():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True)
    parser.add_argument("-o", "--output", type=str, required=True)
    return parser.parse_args()


def main():
    args = get_args()
    output = pl.Path(args.output)
    _LOGGER.info(f"Making output dir {str(output.parents[0])}")
    output.parents[0].mkdir(exist_ok=True, parents=True)
    _LOGGER.info(f"Loading adata from {args.input}")
    adata = read_anndata(args.input)
    _LOGGER.info(f"Computing PCA")
    sc.pp.highly_variable_genes(adata, n_top_genes=4000, subset=True, batch_key="sample")
    x_pca = sc.tl.pca(adata.X)
    latent = pd.DataFrame(x_pca, index=adata.obs_names)
    latent.to_csv(args.output)


if __name__ == '__main__':
    main()
