# Curated Cancer Cell Atlas (CCCA) Data Download Guide

The CanSig-benchmark project provides functionality to download and process datasets from the Curated Cancer Cell Atlas (CCCA). This guide explains how to download datasets and generate processed samples for analysis.

## Overview

The CCCA fetcher (`ccafetcher.py`) handles:
- Downloading raw data and metadata from Dropbox URLs
- Processing ZIP files containing expression matrices
- Converting data to AnnData format for analysis
- Merging metadata with expression data

## Dataset Metadata Structure

Datasets are defined in CSV files (e.g., `3ca.csv`) with the following columns:

- **Tissue**: Tissue type (e.g., lung, breast, brain)
- **Title**: Study author and year
- **Article**: DOI/URL to the original publication
- **Data**: Dropbox URL for the expression data (ZIP format)
- **Meta_data**: Dropbox URL for the metadata (CSV format)
- **Disease**: Disease type or condition
- **Technology**: Sequencing technology (10x, SmartSeq2, etc.)
- **N_cells**: Number of cells in the dataset
- **N_sampels**: Number of samples in the dataset

## Basic Usage

### 1. Download Only (Fetch Mode)

To download raw datasets without processing:

> ⚠️ **Warning:** This will download all datasets specified in the metadata file. Please subset them to the files needed!

```bash
python ccafetcher.py --metadata 3ca.csv --download-dir data/downloads --fetch
```

This will:
- Create the download directory if it doesn't exist
- Download all ZIP files and metadata CSV files
- Store files with hashed names based on their URLs

### 2. Process Single Dataset

To download and convert a specific dataset (first row in the metadata) to AnnData format:

```bash
python ccafetcher.py --metadata 3ca.csv --download-dir data/downloads --dataset-dir data/raw --sample 0
```

Parameters:
- `--metadata`: Path to the metadata CSV file
- `--download-dir`: Directory to store raw downloaded files
- `--dataset-dir`: Directory to store processed AnnData files
- `--sample`: Index of the dataset to process (0-based)
- `--fetch`: Download files if not already present

## Output Format

Processed datasets are saved as `.h5ad` files (AnnData format) in the specified dataset directory. The naming convention is:

- Single dataset per ZIP: `{Author}_{Tissue}.h5ad`
- Multiple datasets per ZIP: `{Author}_{Tissue}_{Subdirectory}.h5ad`

## Data Structure

Each processed AnnData object contains:

- **X**: Expression matrix (sparse CSR format)
- **obs**: Cell-level metadata (observations)
- **var**: Gene-level metadata (variables)
- **uns**: Unstructured metadata including count type

### Changes made to the 3CA.
We corrected a couple of errors in the 3CA.

- Fix the metadata of Neftel et al. by changing the technology to 10x for some samples.
- Converted the .rds in Couturier et al. to .mtx.
- The dataset from Couturier et al. was subsetted to cells that where identified as IDH WT in
genetic_hormonal_features and in histology as GBM.
- The dataset from Yuan et al. was subsetted to GBM patients.