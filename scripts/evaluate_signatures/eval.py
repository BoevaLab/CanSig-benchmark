import argparse
import logging
import pathlib
from typing import List

import numpy as np
import pandas as pd

logging.basicConfig()
logging.root.setLevel(logging.INFO)

_LOGGER = logging.getLogger(__name__)

_SIG_TYPE = "signature_type"
_META_SIG = "meta_signature"
_GROUND_TRUTH = "ground_truth_signature"

def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", "-i", nargs="+",
                        help="List of path to meta signatures stored in .csv.")
    parser.add_argument("--output", "-o", type=str, help="Final results.")
    args = parser.parse_args()
    return args


def get_score(corr: pd.DataFrame, sig_gts: List[str], sig_names: List[str]) -> float:

    corr = corr.copy()
    corr[corr < 0.] = 0.
    contribution = pd.DataFrame(index=sig_names, columns=sig_gts)
    for sig_name in sig_names:
        for sig_gt in sig_gts:
            other_gts = sig_gts.copy()
            other_gts.remove(sig_gt)
            postive_score = corr.loc[sig_name, sig_gt]
            negative_score = np.max(np.maximum(0, corr.loc[[sig_name], other_gts].values -
                                               corr.loc[[sig_gt], other_gts].values))
            score = postive_score - negative_score
            score = np.maximum(0, score)
            contribution.loc[sig_name, sig_gt] = score
    score = 0
    for _ in sig_gts:
        if contribution.shape[0] == 0:
            break
        max_idx = np.unravel_index(np.argmax(contribution.values, axis=None), contribution.values.shape)
        score += contribution.iloc[max_idx]
        contribution = contribution.drop(index=contribution.index[max_idx[0]], columns=contribution.columns[max_idx[1]])

    return score / len(sig_gts)


def main() -> None:
    args = get_args()
    results = pd.DataFrame(index=["scores"])
    for n_cluster, corr_path in enumerate(args.input):
        corr_path = pathlib.Path(corr_path)
        corr_df = pd.read_csv(corr_path, index_col=0)
        sig_type = corr_df.pop(_SIG_TYPE)
        gt_names = corr_df.index[sig_type==_GROUND_TRUTH].tolist()
        meta_sigs_names = corr_df.index[sig_type==_META_SIG].tolist()


        score = get_score(corr_df, gt_names, meta_sigs_names)
        results[n_cluster] = score
    results.to_csv(args.output)


if __name__ == "__main__":
    main()
