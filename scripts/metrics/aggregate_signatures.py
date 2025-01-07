import argparse
import pathlib

import pandas as pd

_CLUSTER = "Cluster"
_SIGNATURE = "Signature"
_SPLIT = "Split"

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", nargs="+")
    parser.add_argument("--output", "-o", type=str)
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    results = []
    for path in args.input:
        path = pathlib.Path(path)
        result = pd.read_csv(path, index_col=0).transpose()
        result.index.name = _CLUSTER
        result[_SIGNATURE] = path.stem.capitalize()
        result[_SPLIT] = path.parent.stem.capitalize()
        result = result.set_index(_SIGNATURE, append=True)
        result = result.reorder_levels([_SIGNATURE, _CLUSTER])

        results.append(result)

    results = pd.concat(results)
    results.to_csv(args.output)


if __name__ == "__main__":
    main()
