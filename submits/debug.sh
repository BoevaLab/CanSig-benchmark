#!/bin/bash
#SBATCH --job-name=snakemake
#SBATCH --time=12:00:00
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G


srun snakemake  --use-singularity --singularity-args "\\-\\-nv" --profile profiles -s integration.smk --keep-going -j 10 --configfile $@