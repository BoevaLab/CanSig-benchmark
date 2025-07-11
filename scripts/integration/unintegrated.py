import logging
import pathlib as pl

import pandas as pd
import scanpy as sc

from utils import read_adata, get_args, set_seed



logging.basicConfig()
logging.root.setLevel(logging.INFO)

_LOGGER = logging.getLogger(__name__)


_ARGS = [("-n", "--n-hvg", str),
         ("-b", "--batch-key", str)]

def main():
    args = get_args(args=_ARGS)
    set_seed(args.random_seed)
    adata = read_adata(args.input, n_hvg=args.n_hvg, batch_key=args.batch_key)
    _LOGGER.info(f"Computing PCA")
    x_pca = sc.tl.pca(adata.X)
    latent = pd.DataFrame(x_pca, index=adata.obs_names)
    latent.to_csv(args.output)


if __name__ == '__main__':
    main()
