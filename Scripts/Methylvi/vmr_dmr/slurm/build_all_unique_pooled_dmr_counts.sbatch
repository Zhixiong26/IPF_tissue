#!/usr/bin/env bash
#SBATCH --job-name=mvi_pdmr_count
#SBATCH --account=lijia
#SBATCH --partition=cpu
#SBATCH --nodelist=cu03
#SBATCH --cpus-per-task=40
#SBATCH --mem=48G
#SBATCH --time=7-00:00:00
#SBATCH --output=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/vmr_dmr_count_%j.out
#SBATCH --error=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/vmr_dmr_count_%j.err

set -Eeuo pipefail
project=/home/lijia/luozhixiong/IPF_tissue
root="$project/Results/MethylVI_vmr_plus_all_unique_pooled_DMR"
bed="$root/shared/dmr/all_unique_hypo_DMRs.merged.bed"
input="$root/shared/dmr/all_unique_hypo_DMRs.input.h5mu"
work="$root/shared/dmr/count_rows"
allc="$project/Results/MethylVI_clean6264_10k30k/vmr/var_0.01/shared/input/selected_cells.allc.tsv"
annotation="$project/Results/Scanpy/E_CYL_ZCP_notebook/cell_id_cell_type.tsv"
[[ -f "$root/shared/dmr/prepare.COMPLETE" ]] || { echo "DMR preparation incomplete" >&2; exit 3; }
features=$(wc -l < "$bed")
(( features >= 2 )) || { echo "Too few merged DMRs: $features" >&2; exit 4; }
cd "$project"
VMR_FILTERED_BED="$bed" VMR_ALLC_TABLE="$allc" VMR_ANNOTATION="$annotation" \
VMR_COUNT_ROWS="$work" VMR_MVI_INPUT="$input" VMR_THREADS="${SLURM_CPUS_PER_TASK:-40}" \
VMR_MIN_COVERED_PERCENT=0 VMR_TARGET_FEATURES="$features" \
  bash "$project/Scripts/Methylvi/vmr/run.sh" build
touch "$root/shared/dmr/counts.COMPLETE"
