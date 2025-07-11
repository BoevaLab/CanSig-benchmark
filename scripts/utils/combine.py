import argparse
import pathlib

import pandas as pd
import numpy as np

batch_mixing_cols = ["kbet", "kbet_per_label"]


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-cluster", nargs="+")
    parser.add_argument("--input-metrics", nargs="+")
    parser.add_argument("--output", type=str)
    return parser.parse_args()


def main():
    args = get_args()
    clustering_results = []
    metrics_results = []
    for file in args.input_cluster:
        file = pathlib.Path(file)
        idx = tuple(str(file).split("/")[-5:-1])
        clustering_result = pd.read_csv(file, index_col=[0, 1])
        mean_clustering_result = clustering_result.mean(1)
        best_resolution = mean_clustering_result.groupby(level=1).mean().idxmax()
        clustering_best = mean_clustering_result.loc[mean_clustering_result.index.get_level_values(level=1)==best_resolution]
        clustering_best = pd.DataFrame(clustering_best.droplevel(1)).transpose()
        clustering_best.index = pd.MultiIndex.from_tuples([idx])
        clustering_results.append(clustering_best)

    for file in args.input_metrics:
        file = pathlib.Path(file)

        idx = tuple(str(file).split("/")[-5:-1])

        metric_result = pd.read_csv(file, index_col=0)
        metric_result.index = pd.MultiIndex.from_tuples([idx])
        metrics_results.append(metric_result)

    clustering_results = pd.concat(clustering_results, axis=0)
    metric_results = pd.concat(metrics_results, axis=0)

    results = pd.concat((clustering_results, metric_results), axis=1)
    results.to_csv(args.output)

if __name__ == '__main__':
    main()
