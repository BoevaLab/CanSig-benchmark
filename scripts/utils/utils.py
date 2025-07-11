from typing import Optional
import scanpy as sc
from pathlib import Path
from argparse import ArgumentParser
import logging
import pathlib
import numpy as np
import pandas as pd

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

_LOGGER = logging.getLogger(__name__)

_SAMPLE = "sample"


def read_adata(path: Path, n_hvg: Optional[int] = 4000, batch_key: bool = False, remove_cycling: bool = True) -> sc.AnnData:
    _LOGGER.info(f"Loading adata from {path}")
    adata = sc.read_h5ad(path)
    _LOGGER.info(f"Loaded adata with {adata.n_obs} cells.")
    if remove_cycling:
        adata = adata[adata.obs["phase"] == "G1"].copy()
        _LOGGER.info(f"Subsetted to {adata.n_obs} none-cycling cells.")
    if n_hvg:
        adata = adata[:, adata.var[f"highly_variable_{n_hvg}_{batch_key}"]].copy()
        _LOGGER.info(f"Subsetted to {adata.n_vars} genes.")        
    return adata


def get_args(args: Optional[list]=None):
    
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", type=Path, required=True)
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("-r", "--random-seed", type=int, required=True)

    if args is not None:
       for short, key, arg_type in args:
            parser.add_argument(short, key, type=arg_type, required=True)
    
    parsed_args = parser.parse_args()
    _LOGGER.info(f"Arguments parsed: {parsed_args}")
    return parsed_args
    

def set_seed(seed):
    import os
    import numpy as np
    import random
    import importlib.util

    
    os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
    np.random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    # Old # No determinism as nn.Upsample has no deterministic implementation
    if importlib.util.find_spec("torch") is not None:
        import torch
        torch.use_deterministic_algorithms(True)
        torch.manual_seed(seed)
        torch.cuda.manual_seed(seed)
        torch.backends.cudnn.enabled = False
        torch.backends.cudnn.benchmark = False
        torch.backends.cudnn.deterministic = True
    if importlib.util.find_spec("scvi") is not None:
        import scvi
        scvi.settings.seed = seed
    _LOGGER.info(f"Set random seed to {seed}")

def save_latent(path: str, latent: np.ndarray, index: np.array):
    out_path = pathlib.Path(path)
    out_path.parent.mkdir(exist_ok=True, parents=True)
    _LOGGER.info(f"Saving latents to {path}")
    if isinstance(latent, dict):
        np.savez(out_path.with_suffix(''), **latent)
    else:
        pd.DataFrame(latent, index=index).to_csv(out_path)