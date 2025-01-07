import argparse
import logging
import pathlib
from typing import List

import numpy as np
import pandas as pd


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", "-i", nargs="+",
                        help="List of path to meta signatures stored in .csv.")
    parser.add_argument("--output", "-o", type=str, help="Final results.")
    args = parser.parse_args()
    return args

def get_score(overlap: pd.DataFrame) -> float:
    contribution = overlap.copy()

    score = 0
    for _ in overlap.columns:
        if contribution.shape[0] == 0:
            break
        max_idx = np.unravel_index(np.argmax(contribution.values, axis=None), contribution.values.shape)
        score += contribution.iloc[max_idx]
        contribution = contribution.drop(index=contribution.index[max_idx[0]], columns=contribution.columns[max_idx[1]])
    return score / overlap.shape[1]

def main() -> None:
    args = get_args()
    results = pd.DataFrame(index=["scores"])
    for overlap_path in args.input:
        overlap_path = pathlib.Path(overlap_path)
        n_cluster = overlap_path.parent.stem
        overlap = pd.read_csv(overlap_path, index_col=0)

        score = get_score(overlap)
        results[n_cluster] = score
    results.to_csv(args.output)


if __name__ == "__main__":
    main()