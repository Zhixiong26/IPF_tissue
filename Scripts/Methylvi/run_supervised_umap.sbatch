#!/usr/bin/env bash
#SBATCH --job-name=ipf_mvi_supumap
#SBATCH --account=lijia
#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/lijia/luozhixiong/IPF_tissue/supervised_umap_%j.out
#SBATCH --error=/home/lijia/luozhixiong/IPF_tissue/supervised_umap_%j.err

set -euo pipefail
export IPF_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd /home/lijia/luozhixiong/IPF_tissue
bash Scripts/03_run_methylvi.sh supervised
