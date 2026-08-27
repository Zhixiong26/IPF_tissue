#!/usr/bin/env bash
#SBATCH --job-name=ipf_vmr_mvi
#SBATCH --account=lijia
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=120G
#SBATCH --time=7-00:00:00
#SBATCH --output=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/vmr_mvi_%j.out
#SBATCH --error=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/vmr_mvi_%j.err

set -euo pipefail
export VMR_THREADS="$SLURM_CPUS_PER_TASK"
export VMR_NUMERIC_THREADS=1
cd /home/lijia/luozhixiong/IPF_tissue
bash Scripts/Methylvi/vmr/run.sh all
