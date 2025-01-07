# Getting Started

## Downloading the package


CanSig can be downloaded directly from github.

```bash
   git clone XXX
   cd CanSig-benchmark
``` 




## Installing the enviorments

This repository provides conda environments for reproducing our benchmarking. We provide two installation options for different reproducibility needs.

**Option 1:** Exact Reproducibility
Use environments from `envs/with_build/`:

```bash

   conda env create -f envs/with_build/CanSig-R.yml
   conda env create -f envs/with_build/CanSig-python.yml
   conda env create -f envs/with_build/cansig-benchmark.yml
```

**Option 2:** Flexible installation

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

## Running the snakemake pipeline


```bash

   snakemake --configfile <path_to_your_config> -c <number_of_cores>
```
To reproduce the results from the paper run

```bash

   snakemake --configfile config/config_gbm.yml -c <number_of_cores>
   snakemake --configfile config/config_breast.yml -c <number_of_cores>
   snakemake --configfile config/config_scc.yml -c <number_of_cores>
   snakemake --configfile config/config_luad.yml -c <number_of_cores>
```

