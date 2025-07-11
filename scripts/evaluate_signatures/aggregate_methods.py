import pathlib

import pandas as pd
import argparse
import pathlib as pl


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", nargs="+")
    parser.add_argument("--output", "-o", type=str)
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    scores = []
    for path in args.input:
        path = pathlib.Path(path)
        score = pd.read_csv(path, index_col=0)
        score.index = [path.parents[1].stem]
        scores.append(score)

    scores = pd.concat(scores)
    scores.index.name = "Method"

    scores.to_csv(args.output)


if __name__ == "__main__":
    main()