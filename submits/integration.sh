snakemake --sdm apptainer --apptainer-args "\\-\\-nvccli"  --configfile configs/integration/umap.yaml -c 24 \
    --resources gpu=1  -s integration.smk --keep-going --rerun-triggers mtime 
snakemake --sdm apptainer --apptainer-args "\\-\\-nvccli"  --configfile configs/integration/batch_key.yaml -c 24 \
    --resources gpu=1  -s integration.smk --keep-going --rerun-triggers mtime 
snakemake --sdm apptainer --apptainer-args "\\-\\-nvccli"  --configfile configs/integration/n_hvgs.yaml -c 24 \
    --resources gpu=1  -s integration.smk --keep-going --rerun-triggers mtime 
snakemake --sdm apptainer --apptainer-args "\\-\\-nvccli"  --configfile configs/integration/umap_wu.yaml -c 24 \
    --resources gpu=1  -s integration.smk --keep-going --rerun-triggers mtime 