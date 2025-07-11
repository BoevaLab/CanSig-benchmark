import pathlib as pl
from argparse import ArgumentParser
from typing import Literal, Tuple

import numpy as np
import pandas as pd

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
        if len(intersection) == 0:
            raise ValueError(
                "There is no common gene between the metasignature and the bulk DataFrame"
            )
        else:
            df = bulk_values.loc[:, intersection]
            if self.std:
                df = (df - df.mean()) / df.std()
            return df.mean(axis=1)


def get_data(
        bulk_file: pl.Path,
        metasignature_file: pl.Path,
        truesignature_file: pl.Path,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Helper function to download files

    Args:

        bulk_file: path to file with the bulk gex
        purity_file: path to the file with the purity info
        metasignature_file: path to the file with the metasignature genes
        truesignature_file: path to the file with the true signature genes

    Returns:

        a scorer object

    """

    bulk_values = pd.read_csv(bulk_file, index_col=0)
    bulk_values = bulk_values.loc[~bulk_values.index.duplicated()]

    metasignatures = pd.read_csv(metasignature_file)
    truesignatures = pd.read_csv(truesignature_file)

    return bulk_values, metasignatures, truesignatures


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
        bulk_file: pl.Path,
        metasignature_file: pl.Path,
        truesignature_file: pl.Path
) -> pd.DataFrame:
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

    bulk_values, metasignatures, truesignatures = get_data(
        bulk_file=bulk_file,
        metasignature_file=metasignature_file,
        truesignature_file=truesignature_file,
    )

    all_scores = get_all_scores(
        bulk_values=bulk_values,
        metasignatures=metasignatures,
        truesignatures=truesignatures,
    )
    return all_scores


def get_args():
    parser = ArgumentParser()
    parser.add_argument("--data-path", "-d", type=str)
    parser.add_argument("--annotation-path", "-a", type=str)
    parser.add_argument("--metasig-path", "-m", type=str)

    return parser.parse_args()


def main():
    args = get_args()
    scores = score_dataset(bulk_file=args.data_path,
                           truesignature_file=args.annotation_path,
                           metasignature_file=args.metasig_path)


if __name__ == '__main__':
    main()
