import scanpy as sc

from utils import get_args, save_latent, read_adata, set_seed
import logging

_SAMPLE_KEY = "sample"
_LATENT = "latent"

logging.basicConfig()
logging.root.setLevel(logging.INFO)

_LOGGER = logging.getLogger(__name__)


_ARGS = [("-n", "--n-hvg", str),
         ("-b", "--batch-key", str)]

def main():
    args = get_args(args=_ARGS)
    set_seed(args.random_seed)
    adata = read_adata(args.input, n_hvg=args.n_hvg, batch_key=args.batch_key)
    sc.tl.pca(adata)
    sc.pp.combat(adata, _SAMPLE_KEY)
    save_latent(args.output, adata.X, adata.obs_names)


if __name__ == '__main__':
    main()
