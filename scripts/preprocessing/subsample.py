
import logging
from argparse import ArgumentParser, Namespace
import scanpy as sc
import numpy as np

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

    if args.n_samples == _ALL:
        _LOGGER.info(f"Didn't subsample the data and wrote the full data to {args.output}.")
        adata.write_h5ad(args.output)
        return

    n_samples = int(args.n_samples)
    assert n_samples <= adata.obs[_SAMPLE].nunique(), "Too few samples in the dataset"
    
    samples = adata.obs[_SAMPLE].unique().tolist()
    _LOGGER.info(f"Initialized generator with random seed {args.random_seed}.")
    rng = np.random.default_rng(seed=args.random_seed)
    samples_to_keep = rng.choice(samples, size=n_samples, replace=False)
    
    bdata = adata[adata.obs[_SAMPLE].isin(samples_to_keep)].copy()
    assert n_samples == bdata.obs[_SAMPLE].nunique(), "Something went wrong, too few samples sampled."
    
    bdata.write_h5ad(args.output)

if __name__ == '__main__':
    main()