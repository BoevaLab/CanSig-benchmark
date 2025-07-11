import logging
import pathlib as pl
import anndata
import pandas as pd
import scanpy as sc
import snapatac2 as snap
from utils import get_args, set_seed, read_adata

_LATENT = "latent"

logging.basicConfig()
logging.root.setLevel(logging.INFO)

_LOGGER = logging.getLogger(__name__)


def get_latent(adata: anndata.AnnData) -> pd.DataFrame:
    BATCH_KEY = "sample"
    snap.tl.spectral(adata, features=None)
    snap.pp.mnc_correct(adata, batch=BATCH_KEY, key_added=_LATENT)
    latent = pd.DataFrame(adata.obsm[_LATENT], index=adata.obs_names)
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
