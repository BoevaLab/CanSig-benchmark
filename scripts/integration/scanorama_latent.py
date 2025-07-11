import numpy as np
import scanpy as sc

from utils import read_adata, get_args, set_seed, save_latent
import logging


logging.basicConfig()
logging.root.setLevel(logging.INFO)

_LOGGER = logging.getLogger(__name__)

_SAMPLE_KEY = "sample"
_LATENT = "latent"

_ARGS = [("-n", "--n-hvg", str),
         ("-b", "--batch-key", str)]

def main():
    args = get_args(args=_ARGS)
    set_seed(args.random_seed)
    adata = read_adata(args.input, n_hvg=args.n_hvg, batch_key=args.batch_key)
    _LOGGER.info("Computing PCA.")
    sc.tl.pca(adata)
    idx = adata.obs_names.copy()
    sorted_idx = np.argsort(adata.obs[_SAMPLE_KEY])
    sorted_adata = adata[sorted_idx].copy()
    sc.external.pp.scanorama_integrate(sorted_adata, _SAMPLE_KEY, adjusted_basis=_LATENT)
    adata = sorted_adata[idx].copy()
    save_latent(args.output, adata.obsm[_LATENT], adata.obs_names)

if __name__ == '__main__':
    main()
