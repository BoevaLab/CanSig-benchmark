import dataclasses

import anndata as ad
import numpy as np
import pytorch_lightning as pl
import torch
import torch.utils.data

import dhaka.vae as vae
import dhaka.data as data
from utils import read_adata, get_args, set_seed, save_latent


_DEFAULT_KEY_ADDED: str = "X_dhaka"


def tensor_to_numpy(t: torch.Tensor) -> np.ndarray:
    return t.detach().cpu().numpy()


_LATENT = "latent"

_ARGS = [("-n", "--n-hvg", str),
         ("-b", "--batch-key", str)]


def run_dhaka(adata: ad.AnnData, key_added: str = _DEFAULT_KEY_ADDED) -> ad.AnnData:
    """Runs the Dhaka algorithm on the specified dataset, learning the representations.

    Args:
        adata: raw gene expression data, will be modified in-place
        config: config
        key_added: field name with representations (to be added to `adata`)

    Returns:
        `adata` with representations added to `obsm[key_added]`
    """

    # Construct the dataset and a training dataloader
    dataset = data.NumpyArrayDataset(
       adata.X.toarray()
    )
    train_dataloader = torch.utils.data.DataLoader(dataset=dataset, shuffle=True, batch_size=64)

    # Define a model and the training procedure (with gradient clipping)
    model = vae.Dhaka(
        n_genes=dataset.n_features,
        latent_dim=3,
        learning_rate=1e-4,
        scale_reconstruction_loss=True,
    )
    trainer = pl.Trainer(
        max_epochs=5,
        gradient_clip_val=2.0,
        gradient_clip_algorithm="norm"
    )
    trainer.fit(model, train_dataloaders=train_dataloader)

    # Get the predictions
    prediction_dataloader = torch.utils.data.DataLoader(dataset=dataset, shuffle=False, batch_size=64)
    model.eval()

    representations_all = np.concatenate([
        tensor_to_numpy(model.representations(batch)) for batch in prediction_dataloader
    ])
    assert representations_all.shape == (adata.shape[0], 3)

    # Add the representations to the AnnData object and return it
    adata.obsm[key_added] = representations_all
    return adata




def main():
    args = get_args(args=_ARGS)
    set_seed(args.random_seed)
    adata = read_adata(args.input, n_hvg=args.n_hvg, batch_key=args.batch_key)
    adata = run_dhaka(adata)
    save_latent(args.output, adata.obsm[_DEFAULT_KEY_ADDED], adata.obs_names)



if __name__ == "__main__":
    main()