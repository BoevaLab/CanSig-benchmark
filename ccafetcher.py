import hashlib
import logging
import pathlib as pl
import warnings
import zipfile
from argparse import ArgumentParser
from itertools import product
from typing import Union

import anndata as anndata
import numpy as np
import pandas as pd
import requests
from scipy.io import mmread
from scipy.sparse import csr_matrix
from tqdm import tqdm

logging.basicConfig()
logging.root.setLevel(logging.INFO)
_LOGGER = logging.getLogger(__name__)

Pathlike = Union[str, pl.Path]

SAMPLE = "sample"


def download(url: str, fname: str, chunk_size=1024):
    resp = requests.get(url, stream=True)
    if resp:
        total = int(resp.headers.get('content-length', 0))
        with open(fname, 'wb') as file, tqdm(
                desc=fname,
                total=total,
                unit='iB',
                unit_scale=True,
                unit_divisor=1024,
        ) as bar:
            for data in resp.iter_content(chunk_size=chunk_size):
                size = file.write(data)
                bar.update(size)
    else:
        raise ValueError(f"Download from {url} failed")

class AnndataDir:
    def __init__(self, path: Pathlike):
        self.path = path
        self.obs_path = self.get_obs_path()
        self.var_path = self.get_var_path()
        self.counts_path = self.get_counts_path()

    @property
    def counts_type(self):
        return self.counts_path.name.split(".")[0]

    def is_valid(self):
        return all([self.obs_path, self.var_path, self.counts_path])

    def get_obs_path(self):
        for file in self.path.iterdir():
            if file.name.endswith(".csv") and "cell" in file.name.lower():
                return file

    def get_var_path(self):
        for file in self.path.iterdir():
            if file.name.endswith(".txt") and "gene" in file.name.lower():
                return file

    def get_counts_path(self):
        for file in self.path.iterdir():
            if file.name.endswith((".rds")):
                _LOGGER.info(f"Found .rds in {str(self.path)}.")
            if file.name.endswith(".mtx"):
                return file

    @staticmethod
    def _read_counts(path) -> csr_matrix:
        if path.name.endswith(".mtx"):
            matrix = mmread(path).transpose()
            matrix = matrix.tocsr()
            matrix = matrix.astype(np.float32)
            return matrix
        else:
            raise NotImplementedError(f"Unknown counts type {path.name}")

    @staticmethod
    def _cast_object_columns_to_string(df: pd.DataFrame) -> pd.DataFrame:
        object_columns = df.select_dtypes(include='object').columns
        if len(object_columns) > 0:
            _LOGGER.info(f"Casting {object_columns.to_list()} to string.")
            df[object_columns] = df[object_columns].fillna("").astype(str)
        return df

    def create_anndata(self, metadata: pd.DataFrame) -> anndata.AnnData:
        if not self.is_valid():
            raise ValueError("Not a valid AnnDataDir.")

        _LOGGER.info("Reading obs.")
        obs = pd.read_csv(self.obs_path.open())
        obs = self.merge_obs_metadata(obs, metadata)
        obs = obs.set_index("cell_name")
        if obs.index.astype(str).str.isdigit().any():
            obs.index = obs.index.astype(str) + "_" + obs[SAMPLE].astype(str)
        obs = self._cast_object_columns_to_string(obs)
        _LOGGER.info("Reading var.")
        var = pd.read_csv(self.var_path.open(), header=None, index_col=0)
        var.index.name = "genes"
        var = self._cast_object_columns_to_string(var)
        _LOGGER.info("Reading counts. This might take a while.")
        counts = self._read_counts(self.counts_path.open())
        counts_type = self.counts_type
        uns = {"count_type": counts_type}
        _LOGGER.info("Creating Anndata.")
        adata = anndata.AnnData(counts, obs=obs, var=var, uns=uns)
        adata.obs_names_make_unique()
        adata.var_names_make_unique()
        return adata

    @staticmethod
    def merge_obs_metadata(obs: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
        metadata = metadata.drop_duplicates()
        if SAMPLE in metadata.columns:
            if obs.index.isin(metadata[SAMPLE]).all():
                _LOGGER.info("Merged the metadata on the cell level.")
                return obs.merge(metadata, how="left", left_index=True,
                                 right_on=SAMPLE, validate="m:1")

        for key_obs, key_metadata in product([SAMPLE, "patient"], [SAMPLE, "patient"]):
            if key_obs in obs.columns and key_metadata in metadata.columns:
                if obs[key_obs].isin(metadata[key_metadata]).all():
                    _LOGGER.info(
                        f"Merging obs and metadata on {key_obs} and {key_metadata}."
                    )
                    return obs.merge(metadata, how="left", left_on=key_obs,
                                     right_on=key_metadata,
                                     suffixes=(None, "_metadata"), validate="m:1")

        raise NotImplementedError("Did not find a way to merge the metadata and obs.")


class CCCAFetcher:
    def __init__(self, meta_data: pd.DataFrame, download_dir: Pathlike,
                 data_column: str = "Data", meta_data_column: str = "Meta_data"):
        self.meta_data_column = meta_data_column
        self.data_column = data_column
        self.meta_data = meta_data
        self.download_dir = pl.Path(download_dir)
        self._make_download_dir()

    def _make_download_dir(self) -> None:
        if self.download_dir.is_dir():
            warnings.warn(f"Download dir {self.download_dir} already exists.")
        self.download_dir.mkdir(parents=True, exist_ok=True)

    def fetch_datasets(self, force_download: bool = False) -> None:
        _LOGGER.info("Start downloading datasets.")
        for name, url, metadata_url in self.meta_data[
            [self.data_column, self.meta_data_column]].itertuples():
            _LOGGER.info(f"Downloading: {name}")
            self.fetch_data(url, force_download=force_download, file_extension="zip")
            self.fetch_data(metadata_url, force_download=force_download,
                            file_extension="csv")

    def fetch_data(self, url: str, force_download: bool, file_extension: str) -> None:
        data_path = self._url_to_data_path(url, file_extension)
        if data_path.is_file() and not force_download:
            return
        download(url, str(data_path))

    def _url_to_data_path(self, url: str, file_extension: str) -> pl.Path:
        url_hash = self._hash_url(url)
        data_path = self.download_dir.joinpath(f"{url_hash}.{file_extension}")
        return data_path

    @staticmethod
    def _hash_url(url: str) -> str:
        hashed_url = hashlib.md5(url.encode("UTF-8"))
        return str(hashed_url.hexdigest())

    def make_anndata(self, data_dir: Pathlike, index: int) -> None:
        data_dir = pl.Path(data_dir)
        data_dir.mkdir(exist_ok=True, parents=True)
        row = self.meta_data.iloc[index]
        name = row.name
        url, metadata_url = row[[self.data_column, self.meta_data_column]]
        _LOGGER.info(f"Start generating AnnData for {name}")
        metadata_path = self._url_to_data_path(metadata_url, file_extension="csv")
        zpath = self._url_to_data_path(url, file_extension="zip")
        if not (metadata_path.is_file() and zpath.is_file()):
            raise ValueError(
                f"Dataset {name} has not been downloaded run fetch_datasets.")
        _LOGGER.info(f"Reading metadata from {metadata_path}.")
        metadata = pd.read_csv(metadata_path)
        _LOGGER.info(f"Reading zipfile from {zpath}.")

        with zipfile.ZipFile(zpath) as zfile:
            root = zipfile.Path(zfile)
            adatadir = AnndataDir(root)
            if adatadir.is_valid():
                _LOGGER.info("Found data at the root level.")
                adata = adatadir.create_anndata(metadata)
                data_path = data_dir.joinpath(f"{name}.h5ad")
                if data_path.is_file():
                    _LOGGER.info(f"No data already exists.")
                    return
                adata.write_h5ad(data_path)

            else:
                _LOGGER.info(
                    "Did not find data at the root level. Checking subdirs.")
                for subdir in root.iterdir():
                    if not subdir.is_dir():
                        continue
                    adatadir = AnndataDir(subdir)
                    data_path = data_dir.joinpath(f"{name}_{subdir.name}.h5ad")
                    if data_path.is_file():
                        _LOGGER.info(f"No data already exists.")
                        continue

                    if not adatadir.is_valid():
                        _LOGGER.info(f"No data found in {subdir.name}")
                        continue
                    _LOGGER.info(f"Found data in a {subdir.name}.")
                    adata = adatadir.create_anndata(metadata)
                    adata.write_h5ad(data_path)


def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--metadata", type=str, required=True)
    parser.add_argument("--download-dir", type=str, required=True)
    parser.add_argument("--dataset-dir", type=str)
    parser.add_argument("--fetch", default=False)
    parser.add_argument("--sample", type=int)
    return parser.parse_args()


# def main():
#     args = parse_args()
#     metadata = pd.read_csv(args.metadata)
#     metadata.index = metadata["Title"].str.split(" ").str[0] + "_" + metadata["Tissue"]
#     fetcher = CCCAFetcher(metadata, args.download_dir)
#     if args.fetch:
#         fetcher.fetch_datasets()
#     if args.sample is not None and args.dataset_dir:
#         fetcher.make_anndata(args.dataset_dir, args.sample)

def main():
    metadata = pd.read_csv("3ca.csv", index_col=0)
    metadata.index = metadata["Title"].str.split(" ").str[0] + "_" + metadata["Tissue"]
    fetcher = CCCAFetcher(metadata, args.download_dir)
    if args.fetch:
        fetcher.fetch_datasets()
    if args.sample is not None and args.dataset_dir:
        fetcher.make_anndata(args.dataset_dir, args.sample)

# def main():
#     metadata = pd.read_csv("3ca.csv", index_col=0)
#     metadata.index = metadata["Title"].str.split(" ").str[0] + "_" + metadata["Tissue"]
#
#     lung_datasets = [f"{author}_lung" for author in ["Bischoff", "Chan", "Kim", "Laughney"]]
#     breast_datasets = [f"{author}_breast" for author in ["Wu", "Pal", "Chung"]]
#     gbm_datasets = [f"{author}_brain" for author in ["Neftel", "Yuan", "Wang"]]
#     other_datasets = ["Ji_skin"]
#     neuroblastoma = [f"{author}_neuroendocrine" for author in ["Dong", "Jansky", "Kildisiute"]]
#     pancrease = [f"{author}_pancreas" for author in ["Hwang", "Peng", "Raghavan", "Steele"]]
#     #all_datasets = lung_datasets + breast_datasets + gbm_datasets + other_datasets
#     all_datasets = breast_datasets[1:2]
#     metadata = metadata.loc[all_datasets]
#     fetcher = CCCAFetcher(metadata, "data/downloads")
#
#     fetcher.fetch_datasets()
#     for i in range(len(all_datasets)):
#         fetcher.make_anndata("data/raw", i)

if __name__ == '__main__':
    main()
