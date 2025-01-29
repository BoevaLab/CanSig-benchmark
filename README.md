# CanSig-benchmark

Complete documentation and tutorials: https://boevalab.github.io/CanSig-benchmark/

## Running the benchmark
### Installing the enviorments
This repository provides conda environments for reproducing our benchmarking. We provide two installation options for different reproducibility needs.

### Option 1: Exact Reproducibility
Use environments from `envs/with_build/`:
```bash
conda env create -f envs/with_build/CanSig-R.yml
conda env create -f envs/with_build/CanSig-python.yml
conda env create -f envs/with_build/cansig-benchmark.yml
```

### Option 2: Flexible Installation
Use environments from `envs/without_build/`:
```bash
conda env create -f envs/without_build/CanSig-R.yml
conda env create -f envs/without_build/CanSig-python.yml
conda env create -f envs/without_build/cansig-benchmark.yml
```

The whole benchmark is run with the cansig-benchmark environment and snakemake activates the other environments. 
```bash
conda activate cansig-benchmark
```

### Downloading the data
To run the benchmark, first download the data from the 3ca.
```
python ccafetcher.py
```
Install the necessary environments.
Then use snakemake to run the benchmark.
```
snakemake --configfile config.yml -c <n-threads> --use-conda
```


### Changes made to the 3CA.
- Fix the metadata of Neftel et al. by changing the technology to 10x for some samples.
- Converted the .rds in Couturier et al. to .mtx.
- The dataset from Couturier et al. was subsetted to cells that where identified as IDH WT in
genetic_hormonal_features and in histology as GBM.
- The dataset from Yuan et al. was subsetted to GBM patients.
