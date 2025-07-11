import logging
import pathlib as pl
import anndata
import pandas as pd
import scanpy as sc

from utils import get_args, set_seed, read_adata

logging.basicConfig()
logging.root.setLevel(logging.INFO)

_LOGGER = logging.getLogger(__name__)


def get_latent(adata: anndata.AnnData) -> pd.DataFrame:
    import scvi
    
    BATCH_KEY = "sample"
    COUNT_LAYER = "counts"

    scvi.model.SCVI.setup_anndata(adata, batch_key=BATCH_KEY, layer=COUNT_LAYER)
    model = scvi.model.SCVI(
        adata)
    _LOGGER.info(f"Start training.")
    model.train(train_size=1.0)
    latent = pd.DataFrame(model.get_latent_representation(), index=adata.obs_names)

    return latent


_ARGS = [("-n", "--n-hvg", str),
         ("-b", "--batch-key", str)]

def main():
    args = get_args(args=_ARGS)
    set_seed(args.random_seed)
    adata = read_adata(args.input, n_hvg=args.n_hvg, batch_key=args.batch_key)
    latent = get_latent(adata)
    _LOGGER.info(f"Writing latent codes to {args.output}")
    latent.to_csv(args.output)


if __name__ == '__main__':
    main()
