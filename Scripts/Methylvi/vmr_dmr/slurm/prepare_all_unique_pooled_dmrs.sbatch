#!/usr/bin/env bash
#SBATCH --job-name=mvi_pdmr_prep
#SBATCH --account=lijia
#SBATCH --partition=cpu
#SBATCH --nodelist=cu03
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/vmr_dmr_prepare_%j.out
#SBATCH --error=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/vmr_dmr_prepare_%j.err

set -Eeuo pipefail
project=/home/lijia/luozhixiong/IPF_tissue
output="$project/Results/MethylVI_vmr_plus_all_unique_pooled_DMR/shared/dmr"
cd "$project"
"/home/lijia/luozhixiong/miniconda3/envs/methylvi/bin/python" \
  Scripts/Methylvi/vmr_dmr/01_prepare_all_unique_pooled_dmrs.py \
  --pairwise-summary Results/Methscan/CYL_ZCP_full_20260826_final/pooled/07_methdiff/pairwise_summary.tsv \
  --blacklist Supplementary/ENCFF356LFX_GRCh38_blacklist.bed.gz \
  --blacklist-fraction 0.2 --raw-p 0.01 --min-abs-diff 0.25 \
  --sort-threads "${SLURM_CPUS_PER_TASK:-8}" --output-dir "$output"
touch "$output/prepare.COMPLETE"
