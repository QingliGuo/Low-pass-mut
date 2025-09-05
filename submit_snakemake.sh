#!/bin/bash
#SBATCH --job-name=lp_mut
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --partition=master-worker
#SBATCH --mem=2000
#SBATCH --time=01:00:00

source ~/.bashrc

cache_path=/data/scratch/DMP/UCEC/GENEVOD/qguo/Software/cache/

export XDG_CACHE_HOME="$cache_path/.cache"
export CONDA_PKGS_DIRS="$cache_path/.conda_pkgs"
export PIP_CACHE_DIR="$cache_path/.pip_cache"
export TMPDIR="$cache_path/.tmp"
export CONDARC="$cache_path/.condarc"

SNAKEMAKE_CONDA_PREFIX="$cache_path/snakemake_conda_envs"

conda activate /data/rds/DMP/UCEC/GENEVOD/qguo/Software/.conda/envs/snakemake_new

snakemake \
  --configfile config.yaml \
  --executor slurm \
  --jobs 50 \
  --use-conda \
  --conda-prefix "$SNAKEMAKE_CONDA_PREFIX" \
  --default-resources slurm_account=dmptxgaag mem_mb=12000 cpus=4 slurm_partition=compute runtime=240 slurm_output="logs/{rule}.{jobid}.out" slurm_error="logs/{rule}.{jobid}.err" \
  --latency-wait 60 \
  --rerun-incomplete \
  --keep-going \
  --printshellcmds \
  --show-failed-logs
