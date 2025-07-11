# Getting started

This repository provides two options for running the Snakemake pipeline: containers hosted on Docker Hub or Conda environments.

## Option 1: Docker Images (Recommended)

### Environment Setup

Install the base environment that includes Snakemake and Apptainer:

```bash
conda env create -f envs/cansig.yml
```
Activate the base environment:

```bash
conda activate cansig
```

### Running the Snakemake Pipeline

Execute the pipeline using Apptainer containers:

```bash
snakemake --configfile <path_to_your_config> --sdm apptainer -c <number_of_cores>
```

### Reproducing Paper Results

To reproduce the results from the paper using Docker containers:

```bash
snakemake --configfile config/config_gbm.yml --sdm apptainer -c <number_of_cores>
snakemake --configfile config/config_breast.yml --sdm apptainer -c <number_of_cores>
snakemake --configfile config/config_scc.yml --sdm apptainer -c <number_of_cores>
snakemake --configfile config/config_luad.yml --sdm apptainer -c <number_of_cores>
```

## Option 2: Conda Environments

### Environment Setup

Use the pre-configured environments from `envs/without_build/`:

```bash
conda env create -f envs/without_build/CanSig-R.yml
conda env create -f envs/without_build/CanSig-python.yml
conda env create -f envs/without_build/cansig-benchmark.yml
```

If you have trouble installing the enviorments with fixed builds, we also provide just the versions at `envs/with_build/`.
The entire benchmark runs within the `cansig-benchmark` environment, while Snakemake automatically activates the other environments as needed:

```bash
conda activate cansig-benchmark
```

### Running the Snakemake Pipeline

Execute the pipeline with Conda environments:

```bash
snakemake --configfile <path_to_your_config> --use-conda -c <number_of_cores>
```

### Reproducing Paper Results

To reproduce the results from the paper, run the following commands:

```bash
snakemake --configfile config/config_gbm.yml -c <number_of_cores>
snakemake --configfile config/config_breast.yml -c <number_of_cores>
snakemake --configfile config/config_scc.yml -c <number_of_cores>
snakemake --configfile config/config_luad.yml -c <number_of_cores>
```

Replace `<number_of_cores>` with the desired number of CPU cores for parallel processing.