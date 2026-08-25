#!/usr/bin/env bash
#SBATCH --job-name=ipf_mvi_umap
#SBATCH --account=lijia
#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/home/lijia/luozhixiong/IPF_tissue/methylvi_umap_%j.out
#SBATCH --error=/home/lijia/luozhixiong/IPF_tissue/methylvi_umap_%j.err

set -euo pipefail
export IPF_THREADS="$SLURM_CPUS_PER_TASK"
export MPLBACKEND=Agg
cd /home/lijia/luozhixiong/IPF_tissue
bash Scripts/03_run_methylvi.sh plots-after
