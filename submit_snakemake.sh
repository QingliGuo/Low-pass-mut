#!/bin/bash
#SBATCH --job-name=snakemake_controller
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1000
#SBATCH --time=02:00:00

source ~/.bashrc
export CONDARC="$PWD/.condarc"
export CONDA_PKGS_DIRS="$PWD/.conda_pkgs"
export PIP_CACHE_DIR="$PWD/.pip_cache"
conda activate /data/rds/DMP/UCEC/GENEVOD/qguo/Software/.conda/envs/snakemake_new

snakemake \
    --configfile config.yaml \
    --executor slurm \
    --jobs 2 \
    --use-conda \
    --conda-frontend mamba \
    --default-resources slurm_account=dmptxgaag mem_mb=16000 cpus=4 slurm_partition=compute runtime=120 jobname={rule}.%j.out \
    --latency-wait 60 \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds \
    --show-failed-logs
