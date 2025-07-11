# Adding your own method

This tutorial explains how to add your method to CanSig.

## Input Data Format

Your method will receive an AnnData object (`.h5ad` file) as input.
- Contains cell phase information in `adata.obs["phase"]`
- Has batch/sample information in `adata.obs["sample"]`
- Raw counts in `adata.X` and normalized, logp1 transformed counts in `adata.layers['normalized']`

## Required Output Format

Your method must output one of two formats:

- A CSV file containing the latent representation
- Each row represents a cell (matching input cells)
- Each column represents a latent dimension
- Row indices should match the cell names from the input AnnData

## Steps to Add Your Method

1. Add your Python Script at `scripts/metasigs/early_integration/your_method_latent.py`:
   ```python
   
   import pathlib as pl
   from argparse import ArgumentParser
   import anndata
   import pandas as pd
   
   def read_anndata(path: str):
       adata = sc.read(path)
       return adata
   
   def get_latent(adata: anndata.AnnData) -> pd.DataFrame:
       # Implement your integration method here
       # Return pd.DataFrame with latent representation

       pass
   
   def get_args():
       parser = ArgumentParser()
       parser.add_argument("-i", "--input", type=str, required=True)
       parser.add_argument("-o", "--output", type=str, required=True)
       # Add any additional arguments your method needs
       return parser.parse_args()
   
   def main():
       args = get_args()
       output = pl.Path(args.output)
       output.parents[0].mkdir(exist_ok=True, parents=True)
       adata = read_anndata(args.input)
       result = get_latent(adata)
       
       # Save based on your method type:
       result.to_csv(output)

   
   if __name__ == '__main__':
       main()
   ```

2. Add Rule to Snakefile:
   ```python
   rule your_method_integration:
       input: rules.subsample.output
       output: cfg.ROOT / "integration/{method}/{scenario}_{n_samples}_{random_seed}/latent.csv"
       benchmark:
           cfg.ROOT / "benchmarks/{scenario}_{n_samples}_{random_seed}/{method}/integration.txt"
       wildcard_constraints: 
           method="your_method.*"  # Adjust pattern as needed
       threads: 4  # Adjust based on your method's requirements
       resources:
           mem_mb=8000,  # Adjust based on your method's requirements
       log:
           cfg.ROOT / "logs/integration/{scenario}_{n_samples}_{random_seed}/{method}/integration.txt"
       shell:
           """python scripts/metasigs/early_integration/your_method_latent.py -i {input} -o {output} {params} &> {log}"""
   ```

3. Add the method to your config file:

```yaml
root: results/your_analysis  # Output directory

methods:
  # List the integration methods you want to use
  your_method:
    param1: value1
    param2: value2
  bbknn:
  harmony:
  genenmf:
  unintegrated:
  scvi:
  scalop:
  cca:

...
```


## Dependencies

If your method requires additional dependencies:
1. Add them to the appropriate environment file
2. Consider creating a separate conda environment if needed
3. Add the `conda:` directive to your Snakefile rule