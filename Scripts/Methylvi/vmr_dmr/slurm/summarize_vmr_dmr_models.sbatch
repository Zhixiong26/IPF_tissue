#!/usr/bin/env bash
#SBATCH --job-name=mvi_vmr_dmr_sum
#SBATCH --account=lijia
#SBATCH --partition=cpu
#SBATCH --nodelist=cu03
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/vmr_dmr_summary_%j.out
#SBATCH --error=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/vmr_dmr_summary_%j.err

set -Eeuo pipefail
project=/home/lijia/luozhixiong/IPF_tissue
root="$project/Results/MethylVI_vmr_plus_all_unique_pooled_DMR"
cd "$project"
/home/lijia/luozhixiong/miniconda3/envs/methylvi/bin/python \
  Scripts/Methylvi/vmr_dmr/03_summarize_vmr_dmr_models.py --root "$root"
touch "$root/workflow.COMPLETE"
