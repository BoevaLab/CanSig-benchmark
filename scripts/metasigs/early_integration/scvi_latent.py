import logging
import pathlib as pl
from argparse import ArgumentParser
import anndata
import pandas as pd
import scanpy as sc

logging.basicConfig()
logging.root.setLevel(logging.INFO)

_LOGGER = logging.getLogger(__name__)



def read_anndata(path: str):
    adata = sc.read(path)
    return adata

def get_latent(adata: anndata.AnnData) -> pd.DataFrame:
    import scvi
    
    BATCH_KEY = "sample"
    COUNT_LAYER = "counts"

    adata = adata[adata.obs["phase"]=="G1"]
    sc.pp.highly_variable_genes(adata, n_top_genes=4000, subset=True, batch_key=BATCH_KEY)

    scvi.model.SCVI.setup_anndata(adata, batch_key=BATCH_KEY, layer=COUNT_LAYER)
    model = scvi.model.SCVI(
        adata)
    _LOGGER.info(f"Start training.")
    model.train(train_size=1.0)
    latent = pd.DataFrame(model.get_latent_representation(), index=adata.obs_names)

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
