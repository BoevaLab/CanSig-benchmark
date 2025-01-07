import logging

import scanpy as sc
import pandas as pd
import argparse
import numpy as np


logging.basicConfig()
logging.root.setLevel(logging.INFO)

_LOGGER = logging.getLogger(__name__)

def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", "-i", type=str,
                        help="Path to meta signatures stored in .csv.")
    parser.add_argument("--output", "-o", type=str, help="Final results.")
    parser.add_argument("--data-path", "-d", type=str, help="Path to the anndata.")
    parser.add_argument("--annotation-path", "-a", type=str,
                        help="Path to the folder holding all the known signatures.")
    args = parser.parse_args()
    return args


def get_overlap(gt_signatures: pd.DataFrame, signatures: pd.DataFrame, var_names: np.array) -> pd.DataFrame:
    results = pd.DataFrame(columns=gt_signatures.columns, index=signatures.columns)

    for gt_sig, gt_list in gt_signatures.items():
        gt_list = set(var_names.intersection(gt_list))
        for name, gene_list in signatures.items():
            gene_list = set(gene_list.dropna())
            results.loc[name, gt_sig] = len(gene_list.intersection(gt_list))

        results[gt_sig] = results[gt_sig] / len(gt_list)
    return results


def main() -> None:
    args = get_args()
    gt_signatures = pd.read_csv(args.annotation_path)
    signatures = pd.read_csv(args.input)
    var_names = sc.read_h5ad(args.data_path).var_names
    overlap = get_overlap(gt_signatures, signatures, var_names)
    overlap.to_csv(args.output)



if __name__ == '__main__':
    main()