#!/usr/bin/env bash
# Build one VMR threshold's shared 30k input and derive its nested 10k H5MU.
# Usage: sbatch build_vmr_feature_inputs.sbatch 0.01|0.02|0.05
#SBATCH --job-name=mvi_vmr_build
#SBATCH --account=lijia
#SBATCH --partition=fat
#SBATCH --nodelist=fat01
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=45
#SBATCH --mem=90G
#SBATCH --time=7-00:00:00
#SBATCH --output=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/vmr_build_%j.out
#SBATCH --error=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/vmr_build_%j.err

set -Eeuo pipefail
if (( $# != 1 )); then echo "Usage: sbatch $0 0.01|0.02|0.05" >&2; exit 2; fi
case "$1" in 0.01|0.02|0.05) ;; *) echo "Unsupported threshold: $1" >&2; exit 3;; esac
project=/home/lijia/luozhixiong/IPF_tissue
source "$project/Scripts/Methylvi/00_experiment_config.sh"
label="var_${1}"
route="$MVI_EXPERIMENT_ROOT/vmr/$label"
shared="$route/shared"
export VMR_METHSCAN_VARIANCE="$1"
export VMR_RESULTS_ROOT="$shared"
export VMR_INPUT_DIR="$shared/input"
export VMR_FILTERED_BED="$shared/input/vmr_filtered.bed"
export VMR_ALLC_TABLE="$shared/input/selected_cells.allc.tsv"
export VMR_COUNT_ROWS="$shared/count_rows"
export VMR_MVI_INPUT="$route/features_30000/input.h5mu"
export VMR_MVI_RESULTS="$route/features_30000/results"
export VMR_TARGET_FEATURES=30000
export VMR_THREADS="${SLURM_CPUS_PER_TASK:-45}"
cd "$project"
bash Scripts/Methylvi/vmr/run.sh prepare
bash Scripts/Methylvi/vmr/run.sh build
"$MVI_METHYLVI_ENV/bin/python" Scripts/Methylvi/shared/06_subset_methylvi_input.py \
  --input "$route/features_30000/input.h5mu" \
  --output "$route/features_10000/input.h5mu" --features 10000
touch "$route/feature_inputs.COMPLETE"
