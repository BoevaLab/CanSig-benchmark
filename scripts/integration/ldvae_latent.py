import scvi

from utils import get_args, save_latent, read_adata, set_seed

_LATENT = "_cansig_latent"

_ARGS = [("-n", "--n-hvg", str),
         ("-b", "--batch-key", str)]

def main():
    args = get_args(args=_ARGS)
    set_seed(args.random_seed)
    adata = read_adata(args.input, n_hvg=args.n_hvg, batch_key=args.batch_key)
    scvi.model.LinearSCVI.setup_anndata(adata, layer="counts", batch_key="sample")
    model = scvi.model.LinearSCVI(adata)
    model.train()
    adata.obsm[_LATENT] = model.get_latent_representation()
    save_latent(args.output, adata.obsm[_LATENT], adata.obs_names)


if __name__ == '__main__':
    main()