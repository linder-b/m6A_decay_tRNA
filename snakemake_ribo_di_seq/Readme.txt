# Snakemake pipeline for analyzing Ribo-seq and disome-seq

# Install enviroment
conda env create -n snakemake_ribo -f environment.yaml
conda activate snakemake_ribo

# adds execute permissions to scripts
chmod +x scripts/*

# prepare STAT and bowtie index
# set up parameters in input.yaml 

# run snakemake pipeline
snakemake --cores 4 --printshellcmds



