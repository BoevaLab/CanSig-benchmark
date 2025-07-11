
import logging
from argparse import ArgumentParser, Namespace
import scanpy as sc
import numpy as np

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

_SAMPLE = "sample"
_ALL = "all"

logging.basicConfig()
logging.root.setLevel(logging.INFO)
_LOGGER = logging.getLogger(__name__)

def get_args() -> Namespace:
    """
    Parse command line arguments for the sample subsetting tool.
    
    Returns:
        Namespace: Parsed command line arguments including:
            - input (str): Path to input h5ad file
            - output (str): Path for output h5ad file
            - n_samples (str): Number of samples to keep or 'all'
            - random_seed (int): Seed for random number generator
    """
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True)
    parser.add_argument("-o", "--output", type=str, required=True)
    parser.add_argument("-n", "--n-samples", type=str, required=True)
    parser.add_argument("-r", "--random-seed", type=int, required=True)
    parser.add_argument("--n-hvgs", nargs='+', type=str, help="Number of highly variable genes to select")
    parser.add_argument("--use-batch-key", nargs="+", help="Batch key to use for HVG selection")
    
    return parser.parse_args()

def main() -> None:
    """
    Main function to execute the sample subsetting workflow.
    
    The function performs the following steps:
    1. Loads an AnnData object from the specified input file
    2. If n_samples is 'all', saves the complete dataset unchanged
    3. Otherwise, randomly selects the specified number of samples
    4. Saves the subsetted data to the specified output file
    
    Raises:
        AssertionError: If requested number of samples exceeds available samples
                       or if the subsetting operation fails to produce expected results
    """
    args = get_args()
    _LOGGER.info(f"Loading data from {args.input}.")
    adata = sc.read_h5ad(args.input)
    samples = adata.obs[_SAMPLE].unique().tolist()

    if args.n_samples != _ALL:
        if float(args.n_samples) <= 1.:
            _LOGGER.info(f"Found {args.n_samples} as percentage")
            n_samples = int(len(samples) * float(args.n_samples))
        else:
            n_samples = int(args.n_samples)
        
        assert n_samples <= len(samples), "Too few samples in the dataset"
        
        _LOGGER.info(f"Initialized generator with random seed {args.random_seed}.")
        rng = np.random.default_rng(seed=args.random_seed)

        
        samples_to_keep = rng.choice(samples, size=n_samples, replace=False)
        
        adata = adata[adata.obs[_SAMPLE].isin(samples_to_keep)].copy()
        _LOGGER.info(f"Subsampled to {adata.obs[_SAMPLE].nunique()} samples.")
        
        assert n_samples == adata.obs[_SAMPLE].nunique(), "Something went wrong, too few samples sampled."
    
    for n_hvg in args.n_hvgs:
        for use_batch_key in args.use_batch_key:
            if n_hvg == "all":
                adata.var[f"highly_variable_{n_hvg}_{use_batch_key}"] = True
            else:
                batch_key = _SAMPLE if use_batch_key == "True" else None
                _LOGGER.info(f"Using {batch_key} to compute hvgs.")
                var = sc.pp.highly_variable_genes(adata, n_top_genes=int(n_hvg), batch_key=batch_key, inplace=False)
                adata.var[f"highly_variable_{n_hvg}_{use_batch_key}"] = var['highly_variable'].values
            _LOGGER.info(f"Added highly_variable_{n_hvg}_{use_batch_key} to adata.var.")

    adata.write_h5ad(args.output)

if __name__ == '__main__':
    main()