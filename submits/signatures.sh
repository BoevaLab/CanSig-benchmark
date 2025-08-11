snakemake --sdm apptainer --apptainer-args "\\-\\-nvccli"  --configfile configs/signatures/gbm.yml -c 8 \
    --resources gpu=1  -s signatures.smk --keep-going --rerun-triggers mtime
snakemake --sdm apptainer --apptainer-args "\\-\\-nvccli"  --configfile configs/signatures/breast.yml -c 8 \
    --resources gpu=1  -s signatures.smk --keep-going --rerun-triggers mtime
snakemake --sdm apptainer --apptainer-args "\\-\\-nvccli"  --configfile configs/signatures/luad.yml -c 8 \
    --resources gpu=1  -s signatures.smk --keep-going --rerun-triggers mtime
snakemake --sdm apptainer --apptainer-args "\\-\\-nvccli"  --configfile configs/signatures/sparse.yml -c 8 \
    --resources gpu=1  -s signatures.smk --keep-going --rerun-triggers mtime
snakemake --sdm apptainer --apptainer-args "\\-\\-nvccli"  --configfile configs/signatures/gbm_subsample.yml -c 8 \
    --resources gpu=1  -s signatures.smk --keep-going --rerun-triggers mtime
snakemake --sdm apptainer --apptainer-args "\\-\\-nvccli"  --configfile configs/signatures/scc.yml -c 8 \
    --resources gpu=1  -s signatures.smk --keep-going --rerun-triggers mtime