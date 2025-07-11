import argparse
import logging
import pathlib as pl
from typing import List, Literal, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

logging.basicConfig()
logging.root.setLevel(logging.INFO)

_LOGGER = logging.getLogger(__name__)

_PHASE = "phase"
_SIG_TYPE = "signature_type"
_META_SIG = "meta_signature"
_GROUND_TRUTH = "ground_truth_signature"
_G1 = "G1"



_NORMMETHOD = Literal["rank", "log_rank", "log", "custom"]
_SCORERNAME = Literal["ssgsea", "average"]


class AvgBulkScorer:
    """This class creates a scorer that takes the average of the (std or not) bulk gex as
    a proxy of signature score"""

    def __init__(self, std: bool):
        """
        Args:

            name: name of the scorer
            std: if True, the bulk gex will be standardized before the average is computed
        Returns:

            None

        """
        self.std = std

    def score(self, bulk_values: pd.DataFrame, metasig: np.ndarray) -> pd.Series:
        """The main scoring function

        Args:

            bulk_values: a df of size (n_samples, n_genes) with the bulk gene expression
            metasig: a list of genes representing the signature to score
        Returns:

            a series with the score for each patient

        """
        intersection = bulk_values.columns.intersection(metasig)
        _LOGGER.info(f"Found {intersection} overlapping genes.")
        if len(intersection) == 0:
            raise ValueError(
                "There is no common gene between the metasignature and the bulk DataFrame"
            )
        else:
            df = bulk_values.loc[:, intersection]
            if self.std:
                df = (df - df.mean()) / df.std()
            return df.mean(axis=1)



def get_all_scores(
        bulk_values: pd.DataFrame,
        metasignatures: pd.DataFrame,
        truesignatures: pd.DataFrame
) -> pd.DataFrame:
    """Computes the scores for all metasignatures and true signatures using
        the appropriate scorer function

    Args:

        bulk_values: df of bulk gex of size (n_samples, n_genes)
        purity: series with the purity information per patient + cancer type of size (n_samples, 2)
        metasignatures: df with n_metasignatures columns, containing in each column the
            list of genes that constitute the metasignature
        truesignatures: df with n_true signatures columns, containing in each column the
            list of genes that constitute the true signature
        std: only used if scorer name is average, if True the bulk gex will be standardized
            before computing the average
        sample_norm_method: only used if scorer name is ssgsea, what method to use
            for sample norm in ssgsea (see ssgsea doc for more info)

    Returns:

        the dataframe of size (n_samples, n_metasignatures + n_true signatures + 1 + 1),
            containing all scores on the metasignatures, the true signatuers, the purity information,
            and TCGA the cancer type the scoring was performed on

    """

    scorer = AvgBulkScorer(std=True)

    all_scores = {}
    for sig in truesignatures.columns:
        all_scores[sig] = scorer.score(
            bulk_values=bulk_values, metasig=truesignatures[sig].ravel()
        )
    for sig in metasignatures.columns:
        all_scores[sig] = scorer.score(
            bulk_values=bulk_values, metasig=metasignatures[sig].ravel()
        )

    all_scores = pd.concat(all_scores, axis=1)

    return all_scores


def score_dataset(
        bulk_data: pd.DataFrame,
        metasignature: pd.DataFrame,
        truesignature: pd.DataFrame
) -> Tuple[List[str], List[str], pd.DataFrame]:
    """Main function, computes the bulk score for metasignatures and reference signatures

    Args:

        bulk_file: path to file with the bulk gex
        metasignature_file: path to the file with the metasignature genes
        truesignature_file: path to the file with the true signature genes
        scorer_name: which scoring to use
        std: only used if scorer name is average, if True the bulk gex will be standardized
            before computing the average
        sample_norm_method: only used if scorer name is ssgsea, what method to use
            for sample norm in ssgsea (see ssgsea doc for more info)

    Returns:

        the dataframe of size (n_samples, n_metasignatures + n_true signatures + 1 + 1),
            containing all scores on the metasignatures, the true signatuers, the purity information,
            and TCGA the cancer type the scoring was performed on

    """

    all_scores = get_all_scores(
        bulk_values=bulk_data,
        metasignatures=metasignature,
        truesignatures=truesignature,
    )
    return truesignature.columns.to_list(), metasignature.columns.to_list(), all_scores



def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", "-i", type=str,
                        help="Path to meta signatures stored in .csv.")
    parser.add_argument("--output", "-o", type=str, help="Final results.")
    parser.add_argument("--data-path", "-d", type=str, help="Path to the anndata.")
    parser.add_argument("--annotation-path", "-a", type=str,
                        help="Path to the folder holding all the known signatures.")
    parser.add_argument("--type", type=str, help="Either scRNA or bulk.")

    args = parser.parse_args()
    return args


def get_adata(data_path: pl.Path) -> ad.AnnData:
    adata = ad.read_h5ad(data_path)
    _LOGGER.info("Subsetting to non-cycling cells.")
    adata = adata[adata.obs[_PHASE] == _G1].copy()

    return adata


def get_scores(adata: ad.AnnData, gt_sigs: pd.DataFrame, meta_sigs: pd.DataFrame):
    gt_names = score_signature(adata, gt_sigs)
    metasig_names = score_signature(adata, meta_sigs)

    return gt_names, metasig_names, adata.obs[gt_names+metasig_names]


def score_signature(adata, sigs: pd.DataFrame) -> List[str]:
    signatures = []
    for sig_name in sigs.columns:
        gene_list = sigs[sig_name].dropna().values[:50]
        if len(adata.var_names.intersection(gene_list))==0:
            continue
        sc.tl.score_genes(adata, gene_list=gene_list, score_name=sig_name)
        signatures.append(sig_name)
    return signatures


def corr_signatures(df: pd.DataFrame, gt_names: List[str],  meta_sigs_names: List[str]) -> pd.DataFrame:
    corr = df[gt_names + meta_sigs_names].corr(method="spearman")
    corr[_SIG_TYPE] = _META_SIG
    corr.loc[corr.index.isin(gt_names), _SIG_TYPE] = _GROUND_TRUTH
    return corr


def main() -> None:
    args = get_args()
    gt_sigs = pd.read_csv(args.annotation_path)
    meta_sigs = pd.read_csv(args.input)

    if args.type.lower() == "scrna":
        adata = get_adata(args.data_path)
        gt_names, meta_sigs_names, scores = get_scores(adata, gt_sigs, meta_sigs)
    elif args.type.lower() == "bulk":
        bulk_data = pd.read_csv(args.data_path, index_col=0)
        gt_names, meta_sigs_names, scores = score_dataset(bulk_data=bulk_data,
                               truesignature=gt_sigs,
                               metasignature=meta_sigs)
    else:
        raise ValueError(f"Unkown scoring type {args.type}.")

    corr = corr_signatures(scores, gt_names, meta_sigs_names)
    corr.to_csv(args.output)



if __name__ == "__main__":
    main()
