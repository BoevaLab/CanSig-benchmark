# Running CanSig on your data

This guide explains how to run CanSig on your own datasets. The pipeline uses a configuration-based approach where you specify your data, methods, and parameters in YAML files.

## Configuration Structure

Add your configuration file to the `configs` folder. It should follow this structure:

```yaml
root: results/your_analysis  # Output directory

methods:
  # List the integration methods you want to use
  bbknn:
  harmony:
  genenmf:
  unintegrated:
  scvi:
  scalop:
  cca:

scenarios:
  your_scenario_name:
    data_path: path/to/your/data.h5ad  # Path to your AnnData file
    signatures:
      - signature_name1
      - signature_name2
    preprocessing:  # Optional preprocessing parameters
      min_genes: 1000  # Minimum number of genes per cell
      max_pct_mt: 20   # Maximum percentage of mitochondrial genes
      excluded_samples: ["sample1", "sample2"]  # Samples to exclude

signatures:
  signature_name1:
    annotation_path: path/to/annotations1.csv
    scoring_scenario: your_scenario_name # Dataset used for comapring gene signature scores.
    n_cluster: 4
 
```

## Input Data
CanSig expects input data in `AnnData` (.h5ad) format with specific requirements for the AnnData object structure. The `adata.obs` must contain a "sample" column that stores the batch ID for each cell in the dataset. The `adata.var_names` should contain gene symbols as the feature identifiers. The `adata.uns` dictionary must include a "counts_type" field specifying the sequencing platform used to generate the data. Supported values for counts_type are "10x" for 10x Genomics data, "microwell array-based platform" for Microwell data, "microwell-seq" for Microwell-seq data, "smartseq2" or "SmartSeq2" for Smart-seq2 data, and "seqwell" for SeqWell data. All of these sequencing technologies are used in the currate cancer cell atlas.


   - Create an entry for each cell type signature
   - Specify the path to annotation files
   - Set the number of expected clusters (`n_cluster`)
   - Link to the appropriate scenario using `scoring_scenario`

## Running the Pipeline

Once your configuration is set up:

1. Place your configuration file in the root directory
2. Run the pipeline using:
   ```bash
   snakemake --configfile your_config.yml -c <number-of-cores>
   ```


## Output Structure

Results will be organized in your specified directory as follows:
```
results/your_analysis/
├── benchmarks/     # Performance benchmarks for each pipeline step
├── corrs/          # Correlation analysis results
├── data/           # Processed data files
├── integration/    # Integrated datasets for each method
├── logs/          # Pipeline execution logs
├── metasigs/      # Generated meta-signatures
└── scores/        # Evaluation scores and metrics
```
Each directory contains the intermediate and final results for that stage of the pipeline.

## Advanced Configuration

### Random Sampling
You can configure random down sampling of samples by adding these parameters to your scenario:

```yaml
scenarios:
  your_scenario:
    n_samples: [5, 10, 15]  # List of sample sizes to try
    random_seeds: 5  # Number of random seeds for replication
```

### Bulk RNA Integration
If you want to comapre your signatures also on bulk data add:

```yaml
scenarios:
  your_scenario:
    bulk_path: path/to/bulk_data.h5ad  # Path to bulk RNA-seq data
```
