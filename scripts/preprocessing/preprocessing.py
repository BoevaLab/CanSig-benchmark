import anndata
import logging
import numpy as np
import pathlib as pl
import scanpy as sc
import warnings
from argparse import ArgumentParser, Namespace
from os import PathLike
from scipy.sparse import csr_matrix, diags
from sklearn.utils import sparsefuncs
from typing import Tuple, Literal, List, Optional
import pandas as pd

warnings.simplefilter(action='ignore', category=FutureWarning)

SEQ_10x = "10x"
MICROWELL = 'microwell array-based platform'
MICROWELL_SEQ = "microwell-seq"
SEQ_SS2 = "smartseq2"
SMARTSEQ2 = "SmartSeq2"
SEQWELL = "seqwell"
SAMPLE = "sample"
NORMALIZED_KEY = "normalized"
COUNT_TYPE = "count_type"

logging.basicConfig()
logging.root.setLevel(logging.INFO)

_LOGGER = logging.getLogger(__name__)


def get_args() -> Namespace:
    """
    Parse command line arguments for the preprocessing pipeline.
    
    Returns:
        Namespace: Parsed command line arguments including:
            - input: Path to input h5ad file
            - output: Path for output h5ad file
            - excluded_sample: List of sample names to exclude
            - min_genes: Minimum number of genes required per cell
            - min_counts: Minimum number of counts required per cell
            - max_pct_mt: Maximum percentage of mitochondrial genes allowed
    """
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True)
    parser.add_argument("-o", "--output", type=str, required=True)
    parser.add_argument("--excluded-sample", type=str, nargs="+")
    parser.add_argument("--min-genes", type=int)
    parser.add_argument("--min-counts", type=int)
    parser.add_argument("--max-pct-mt", type=float, required=True)
    parser.add_argument("--labels-path", type=pl.Path)
    return parser.parse_args()


def get_counts_per_cell(x: csr_matrix):
    """
    Calculate the total count of gene expression for each cell.
    
    Args:
        x (csr_matrix): Sparse matrix containing gene expression data
        
    Returns:
        numpy.ndarray: Array containing total counts for each cell
    """
    return np.asarray(x.sum(1)).ravel()


def get_counts_from_tpm(tpm: csr_matrix, technology: str) -> csr_matrix:
    """
    Convert TPM (Transcripts Per Million) values to estimated count data.
    
    Args:
        tpm (csr_matrix): Sparse matrix of TPM values
        technology (str): Sequencing technology used
        
    Returns:
        csr_matrix: Estimated count data
        
    Raises:
        NotImplementedError: If the sequencing technology is not recognized
    """
    tpm = tpm.copy()
    if any(technology == sequencing_tech for sequencing_tech in [SEQ_10x, MICROWELL, MICROWELL_SEQ, SEQWELL]):
        library_sizes = np.zeros(tpm.shape[0])
        const = get_counts_per_cell(tpm).mean()
        for n_row in range(tpm.shape[0]):
            row = tpm[n_row]
            library_size = const / row[row > 0.0].min()
            library_sizes[n_row] = library_size

        counts = diags(library_sizes / const) * tpm
        counts.data = np.round(counts.data)
        _LOGGER.info(f"Converting TPMs to counts by estimating library size. Estimated library size:\n"
                     f" Mean: {library_sizes.mean()}\n Max: {library_sizes.max()}\n Min: {library_sizes.min()}")
    elif technology == SEQ_SS2:
        _LOGGER.info("Converting TPMs to counts by rounding.")
        counts = tpm.copy()
        counts.data = np.round(counts.data)
    else:
        raise NotImplementedError(f"Unknown technology {technology}.")

    return counts


def normalize(counts: csr_matrix, target_counts=1e5) -> csr_matrix:
    """
    Normalize counts to a target sum per cell.
    
    Args:
        counts (csr_matrix): Raw count data
        target_counts (float): Target sum for each cell after normalization
        
    Returns:
        csr_matrix: Normalized count data
    """
    library_size = get_counts_per_cell(counts)
    sparsefuncs.inplace_row_scale(counts, target_counts / library_size)
    return counts


def get_tpm_counts(input, count_type: str, technology: str) -> Tuple[csr_matrix, csr_matrix]:
    """
    Process input data to get both TPM and count matrices.
    
    Args:
        input: Input count or TPM data
        count_type (str): Type of count data provided
        technology (str): Sequencing technology used
        
    Returns:
        Tuple[csr_matrix, csr_matrix]: TPM values and count data
    """
    _LOGGER.info(f"Found {count_type} count type.")
    if count_type in ['Exp_data_UMIcounts', "Exp_data_UMIcounts_10X"]:
        counts = input.copy()
        tpm = normalize(input)
        return tpm, counts

    counts = get_counts_from_tpm(input, technology)
    tpm = normalize(input)
    return tpm, counts


def set_tpm_counts(adata: anndata.AnnData, technology: str) -> anndata.AnnData:
    """
    Set TPM and count data in AnnData object and perform log transformation.
    
    Args:
        adata (anndata.AnnData): Input data object
        technology (str): Sequencing technology used
        
    Returns:
        anndata.AnnData: Processed data object with TPM and counts
    """
    count_type = adata.uns[COUNT_TYPE]
    tpm, counts = get_tpm_counts(adata.X, count_type, technology)
    del adata.uns[COUNT_TYPE]
    adata.layers["counts"] = counts
    adata.X = tpm
    _LOGGER.info("Log plus one transforming the normalized data (adata.X). ")
    sc.pp.log1p(adata, base=2)
    return adata


def remove_low_count_cells(adata: anndata.AnnData, min_counts: Optional[int],
                           min_genes: Optional[int]) -> anndata.AnnData:
    """
    Filter cells based on minimum count and gene expression thresholds.
    
    Args:
        adata (anndata.AnnData): Input data object
        min_counts (Optional[int]): Minimum counts required per cell
        min_genes (Optional[int]): Minimum genes required per cell
        
    Returns:
        anndata.AnnData: Filtered data object
    """
    if min_counts:
        _LOGGER.info(f"Removing cells with less than {min_counts} counts.")
        sc.pp.filter_cells(adata, min_counts=min_counts)
    if min_genes:
        _LOGGER.info(f"Removing cells with less than {min_genes} genes expressed.")
        sc.pp.filter_cells(adata, min_genes=min_genes)
    _LOGGER.info(f"Kept {adata.n_obs} cells.")
    return adata


def subset_malignant(adata: anndata.AnnData) -> anndata.AnnData:
    """
    Subset data to include only malignant cells.
    
    Args:
        adata (anndata.AnnData): Input data object
        
    Returns:
        anndata.AnnData: Data object containing only malignant cells
        
    Raises:
        ValueError: If malignant cells cannot be identified
    """
    MALIGNANT = "malignant"
    if MALIGNANT in adata.obs.columns:
        _LOGGER.info("Used adata.obs['malignant'] to find malignant cells.")
        idx = adata.obs[MALIGNANT] == "yes"
    elif MALIGNANT in adata.obs["cell_type"].str.lower().values:
        _LOGGER.info(
            f"Used cell_types {adata.obs.cell_type.unique().tolist()} to find malignant cells.")
        idx = adata.obs["cell_type"].str.lower() == MALIGNANT
    else:
        raise ValueError("No way of determining which cells are malignant.")

    adata = adata[idx].copy()
    _LOGGER.info(f"Found {adata.n_obs} malignant cells.")
    return adata


def determine_seq_technology(adata: anndata.AnnData):
    """
    Determine the sequencing technology used from the AnnData object.
    
    Args:
        adata (anndata.AnnData): Input data object
        
    Returns:
        str: Sequencing technology identifier
        
    Raises:
        KeyError: If technology information is missing
        ValueError: If multiple technologies are found
    """
    TECHNOLOGY = "technology"
    if TECHNOLOGY not in adata.obs.columns:
        raise KeyError(f"Key {TECHNOLOGY} not in adata.obs.columns.")
    technology = adata.obs[TECHNOLOGY].unique()
    if len(technology) != 1:
        raise ValueError("More than one technology found.")
    technology = technology[0]
    _LOGGER.info(f"Determined {technology} as sequencing technology.")
    return technology.lower()


def remove_high_mt_cells(adata: anndata.AnnData, max_pct_mt: float) -> anndata.AnnData:
    """
    Remove cells with high mitochondrial gene expression.
    
    Args:
        adata (anndata.AnnData): Input data object
        max_pct_mt (float): Maximum percentage of mitochondrial counts allowed
        
    Returns:
        anndata.AnnData: Filtered data object
    """
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False,
                               inplace=True)
    adata = adata[adata.obs["pct_counts_mt"] <= max_pct_mt].copy()
    _LOGGER.info(f"Kept {adata.n_obs} cells after removing high MT cells.")
    return adata


def remove_samples(adata: anndata.AnnData) -> anndata.AnnData:
    """
    Remove samples with fewer than 50 cells in G1 phase.
    
    Args:
        adata (anndata.AnnData): Input data object
        
    Returns:
        anndata.AnnData: Filtered data object
    """
    obs = adata.obs.copy()
    obs_no_cc = obs[obs["phase"]=="G1"]
    cells_per_sample = obs_no_cc[SAMPLE].value_counts()

    samples_to_keep = cells_per_sample[cells_per_sample > 100].index.tolist()
    adata = adata[adata.obs[SAMPLE].isin(samples_to_keep)].copy()
    _LOGGER.info(f"Removed {cells_per_sample.shape[0] - len(samples_to_keep)} of {cells_per_sample.shape[0]} samples.")
    return adata


def r_names(adata: anndata.AnnData) -> anndata.AnnData:
    """
    Convert variable and observation names to R-compatible format.
    
    Args:
        adata (anndata.AnnData): Input data object
        
    Returns:
        anndata.AnnData: Data object with R-compatible names
    """
    _LOGGER.info("Changing `obs_names` and `var_names` to R names.")
    adata.var_names = adata.var_names.str.replace("_", "-")
    adata.obs_names = adata.obs_names.str.replace("_", "-")
    return adata


def remove_excluded_samples(adata: sc.AnnData, excluded_samples: List[str]) -> sc.AnnData:
    """
    Remove specified samples from the dataset.
    
    Args:
        adata (sc.AnnData): Input data object
        excluded_samples (List[str]): List of sample names to exclude
        
    Returns:
        sc.AnnData: Filtered data object
    """
    _LOGGER.info(f"Anndata contains {adata.obs[SAMPLE].unique().tolist()} as samples.")
    adata = adata[~adata.obs[SAMPLE].isin(excluded_samples)].copy()
    _LOGGER.info(f"After removing excluded samples anndata contains {adata.obs[SAMPLE].unique().tolist()} as samples.")
    return adata


def preprocessing(adata: anndata.AnnData, excluded_samples: List[str], min_genes: int, min_counts: int,
                  max_pct_mt: float) -> anndata.AnnData:
    """
    Perform complete preprocessing pipeline on single-cell RNA sequencing data.
    
    Args:
        adata (anndata.AnnData): Input data object
        excluded_samples (List[str]): List of sample names to exclude
        min_genes (int): Minimum number of genes required per cell
        min_counts (int): Minimum number of counts required per cell
        max_pct_mt (float): Maximum percentage of mitochondrial counts allowed
        
    Returns:
        anndata.AnnData: Fully preprocessed data object
    """
    _LOGGER.info(
        f"Started preprocessing with {adata.n_obs} cells, {adata.n_vars} genes and "
        f"{adata.obs[SAMPLE].nunique()} samples.")

    if "normalized_Exp_data_TPM" == adata.uns[COUNT_TYPE]:
        adata.X.data = (2 ** adata.X.data) - 1

    technology = determine_seq_technology(adata)
    adata = subset_malignant(adata)
    adata = remove_low_count_cells(adata, min_genes=min_genes, min_counts=min_counts)
    adata = remove_high_mt_cells(adata, max_pct_mt)
    adata = score_cell_cycle(adata)
    adata = remove_samples(adata)
    if excluded_samples:
        adata = remove_excluded_samples(adata, excluded_samples)
    adata = set_tpm_counts(adata, technology)
    adata = r_names(adata)
    sc.pp.filter_genes(adata, min_cells=int(0.001 * adata.n_obs))
    _LOGGER.info(
        f"Finished preprocessing with {adata.n_obs} cells, {adata.n_vars} genes and "
        f"{adata.obs[SAMPLE].nunique()} samples.")
    return adata


def read_anndata(path: PathLike) -> anndata.AnnData:
    _LOGGER.info(f"Reading Anndata from {path}.")
    adata = anndata.read_h5ad(path)
    if adata.obs_names.str.isdigit().any():
        _LOGGER.info("Found integers in adata.obs_names appending sample name "
                     "to ensure casting to strings.")
        adata.obs_names = adata.obs_names.astype(str) + "_" + adata.obs[SAMPLE].astype(str)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    return adata


def write_anndata(adata: anndata.AnnData, output_path: PathLike) -> None:
    output_path = pl.Path(output_path)
    output_path.parent.mkdir(exist_ok=True, parents=True)
    adata.strings_to_categoricals()
    _LOGGER.info(f"Writing AnnData to {output_path}.")
    adata.write_h5ad(output_path)


def score_cell_cycle(adata: anndata.AnnData) -> anndata.AnnData:
    s_genes = ['MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG',
               'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP',
               'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76',
               'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2',
               'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2',
               'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8']

    g2m_genes = ['HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A',
                 'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF',
                 'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB',
                 'BUB1', 'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP',
                 'CDCA3', 'HN1', 'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1',
                 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8', 'ECT2', 'KIF23', 'HMMR',
                 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF', 'NEK2',
                 'G2E3', 'GAS2L3', 'CBX5', 'CENPA']

    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    return adata

def add_labels(adata, labels_path) -> sc.AnnData:
    _LOGGER.info(f"Reading labels from {labels_path}")
    labels = pd.read_csv(labels_path, index_col=0)
    _LOGGER.info(f"Adding labels to {adata.n_obs} cells.")
    adata.obs = adata.obs.merge(labels, how="left", left_index=True, right_index=True)
    adata = adata[~adata.obs["celltype"].isna()].copy()
    _LOGGER.info("Added labels to {adata.n_obs} cells.")
    return adata
    


def make_datadir(output):
    output_path = pl.Path(output)
    output_path.parents[0].mkdir(exist_ok=True, parents=True)


def main():
    args = get_args()
    make_datadir(args.output)
    adata = read_anndata(args.input)
    print(args.labels_path)
    if args.labels_path.name.endswith(".csv"):
        adata = add_labels(adata, args.labels_path)
        
    adata = preprocessing(adata, excluded_samples=args.excluded_sample, min_genes=args.min_genes,
                          max_pct_mt=args.max_pct_mt, min_counts=args.min_counts)
    write_anndata(adata, args.output)


if __name__ == '__main__':
    main()