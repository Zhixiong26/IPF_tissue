#!/usr/bin/env bash
#SBATCH --job-name=ipf_allc_smoke
#SBATCH --account=lijia
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/allc_smoke_%j.out
#SBATCH --error=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/allc_smoke_%j.err

set -euo pipefail
export IPF_THREADS="$SLURM_CPUS_PER_TASK"
export IPF_NUMERIC_THREADS=1

cd /home/lijia/luozhixiong/IPF_tissue
bash Scripts/Methylvi/allcools/run.sh smoke
