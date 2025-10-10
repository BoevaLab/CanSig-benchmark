# Reproducing Our Results

This tutorial provides step-by-step instructions for reproducing the benchmark results from the CanSig-benchmark paper using the curated cancer datasets from the Curated Cancer Cell Atlas (CCCA).

## Prerequisites

### Environment Setup
Ensure you have the CanSig environment set up and activated:

```bash
conda env create -f envs/cansig.yml
conda activate cansig
```

### Data Preparation
Before running the pipeline, you need to download the cancer datasets from the CCCA. This step is crucial for reproducing the paper results.

For detailed instructions, see the [CCCA Data Download Guide](../ccca.md).

Quick setup for a subset of datasets:
```bash
python ccafetcher.py --metadata 3ca_test.csv --download-dir data/downloads --dataset-dir data/raw --sample 0 --fetch
```

## Running the Complete Benchmark

### Signature Analysis Pipeline

The signature analysis pipeline evaluates different gene signature methods across cancer datasets:

```bash
snakemake --configfile config/signatures/<config>.yaml --sdm apptainer --apptainer-args "\\-\\-nvccli" -c <number_of_cores> -s signatures.smk
```

Available signature configuration files can be found at `configs/signatures`.

### Integration Analysis Pipeline

The integration analysis pipeline benchmarks batch integration methods:

```bash
snakemake --configfile config/integration/<config>.yaml --sdm apptainer --apptainer-args "\\-\\-nvccli" -c <number_of_cores> -s integration.smk
```

Available signature configuration files can be found at `configs/integration`.

### GPU Support

The `--apptainer-args "\\-\\-nvccli"` argument enables NVIDIA container support for GPU-accelerated methods. This is optional but recommended if you have compatible GPU hardware.

For CPU-only execution, remove this argument:
```bash
snakemake --configfile config/signatures/<config>.yaml --sdm apptainer -c <number_of_cores> -s signatures.smk
```

## Reproducing Specific Results

### Paper Figure Reproduction

After running the pipelines, you can reproduce the figures from the paper using the provided Jupyter notebooks in the `figures/` directory:

```bash
# Navigate to figures directory
cd figures/

# Run figure generation notebooks
jupyter notebook fig2.ipynb    # Main benchmark results
jupyter notebook fig3.ipynb    # Method comparison
jupyter notebook fig4_1.ipynb  # Performance analysis part 1
jupyter notebook fig4_2.ipynb  # Performance analysis part 2
```

### Supplementary Figures

```bash
jupyter notebook supfig3.ipynb  # Supplementary figure 3
jupyter notebook supfig4.ipynb  # Supplementary figure 4
```

## Expected Output Structure

The pipeline will generate results in the following structure:

```
results/
├── signatures/
│   ├── benchmarks/     # Performance benchmarks
│   ├── data/          # Processed signature data
│   ├── logs/          # Execution logs
│   └── scores/        # Signature evaluation scores
└── integration/
    ├── benchmarks/     # Integration benchmarks
    ├── integration/    # Integrated datasets
    ├── logs/           # Execution logs
    └── corrs/          # Correlation analyses
```

