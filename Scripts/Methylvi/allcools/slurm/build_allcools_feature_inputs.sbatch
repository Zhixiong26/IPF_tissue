#!/usr/bin/env bash
# Build the shared 30k ALLCools input and derive its nested 10k H5MU.
#SBATCH --job-name=mvi_allc_build
#SBATCH --account=lijia
#SBATCH --partition=fat
#SBATCH --nodelist=fat01
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=45
#SBATCH --mem=90G
#SBATCH --time=7-00:00:00
#SBATCH --output=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/allc_build_%j.out
#SBATCH --error=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/allc_build_%j.err

set -Eeuo pipefail
project=/home/lijia/luozhixiong/IPF_tissue
source "$project/Scripts/Methylvi/00_experiment_config.sh"
route="$MVI_EXPERIMENT_ROOT/allcools"
shared="$route/shared"
export IPF_MVI_ROOT="$shared"
export IPF_ALLCOOLS_ROOT="$shared/allcools_5kb"
export IPF_ALLC_DIR="$shared/input_allc"
export IPF_ALLC_TABLE="$shared/selected_cells.allc.tsv"
export IPF_MCDS="$shared/mcg_5kb.mcds"
export IPF_ALLCOOLS_H5AD="$route/features_30000/allcools_30k.h5ad"
export IPF_MVI_INPUT="$route/features_30000/input.h5mu"
export IPF_MVI_RESULTS="$route/features_30000/results"
export IPF_TARGET_FEATURES=30000
export IPF_THREADS="${SLURM_CPUS_PER_TASK:-45}"
cd "$project"
bash Scripts/Methylvi/allcools/run.sh prepare
bash Scripts/Methylvi/allcools/run.sh cluster
bash Scripts/Methylvi/allcools/run.sh build
"$MVI_METHYLVI_ENV/bin/python" Scripts/Methylvi/shared/06_subset_methylvi_input.py \
  --input "$route/features_30000/input.h5mu" \
  --output "$route/features_10000/input.h5mu" --features 10000
touch "$route/feature_inputs.COMPLETE"
