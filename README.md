# Getting Started

The CanSig-benchmark is a comprehensive benchmarking framework for cancer gene signature methods. This guide will walk you through setting up the environment and preparing the data needed to run the benchmarking pipeline.

### Environment Setup

The base environment contains Snakemake, Apptainer, and other essential tools needed to orchestrate the benchmarking pipeline. All computational methods are containerized to ensure reproducibility across different systems.

Install the base environment that includes Snakemake and Apptainer:

```bash
conda env create -f envs/cansig.yml
```
Activate the base environment:

```bash
conda activate cansig
```

### Data Preparation

After setting up the environment, you'll need to download the cancer datasets from the Curated Cancer Cell Atlas (CCCA). These datasets form the foundation for the benchmarking analysis.
For detailed instructions on downloading and processing the CCCA datasets, see the [CCCA Data Download Guide](https://boevalab.github.io/CanSig-benchmark/ccca.html).

Quick start for downloading a test dataset:

```bash
python ccafetcher.py --metadata 3ca_test.csv --download-dir data/downloads --dataset-dir data/raw --sample 0 --fetch
```

### Running the Snakemake Pipeline

Once you have set up the environment and downloaded the necessary datasets, you can execute the benchmarking pipeline. The pipeline uses docker containers to ensure all methods run in isolated, reproducible environments.

Execute the pipeline using docker containers:

```bash
snakemake --configfile <path_to_your_config> --sdm apptainer -c <number_of_cores> -s <snakemake_file>
```

Replace `<path_to_your_config>` with your specific configuration file, `<number_of_cores>` with the number of CPU cores you want to allocate and  `<snakemake_file>` with either `signatures.smk` or `integration.smk` depending on which part of the pipeline you wish to run.

### Reproducing Paper Results

To reproduce the specific results from the CanSig-benchmark paper, use the provided configuration files. The pipeline includes two main components: signature analysis and integration analysis.

**Signature Analysis Pipeline:**
```bash
snakemake --configfile config/signatures/<config>.yaml --sdm apptainer --apptainer-args "\\-\\-nvccli" -c <number_of_cores> -s signatures.smk
```

**Integration Analysis Pipeline:**
```bash
snakemake --configfile config/integration/<config>.yaml --sdm apptainer --apptainer-args "\\-\\-nvccli" -c <number_of_cores> -s integration.smk
```

Replace `<config>` with the specific configuration file name you want to use. The `--nvccli` argument enables NVIDIA container support for GPU-accelerated methods.

## Next Steps

1. **Configure your analysis**: Examine the configuration files in `config/signatures/` and `config/integration/` to understand the available options
3. **Analyze results**: Once done running your benchmark you can analyse the created data using the notebooks provided in ```figures/```
