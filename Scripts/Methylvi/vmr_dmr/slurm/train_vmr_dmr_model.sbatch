#!/usr/bin/env bash
# Usage: sbatch train_vmr_dmr_model.sbatch THRESHOLD 10000|30000
#SBATCH --job-name=mvi_vmr_dmr_train
#SBATCH --account=lijia
#SBATCH --partition=cpu
#SBATCH --nodelist=cu03
#SBATCH --cpus-per-task=32
#SBATCH --mem=48G
#SBATCH --time=7-00:00:00
#SBATCH --output=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/vmr_dmr_train_%j.out
#SBATCH --error=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/vmr_dmr_train_%j.err

set -Eeuo pipefail
if (( $# != 2 )); then echo "Usage: $0 0.01|0.02|0.05 10000|30000" >&2; exit 2; fi
case "$1" in 0.01|0.02|0.05) ;; *) exit 3;; esac
case "$2" in 10000|30000) ;; *) exit 4;; esac
project=/home/lijia/luozhixiong/IPF_tissue
route="$project/Results/MethylVI_vmr_plus_all_unique_pooled_DMR/vmr_dmr/var_${1}/features_${2}"
input="$route/input.h5mu"; results="$route/results"
[[ -f "$route/input.COMPLETE" && -s "$input" ]] || exit 5
export IPF_MVI_INPUT="$input" IPF_MVI_RESULTS="$results" VMR_MVI_RESULTS="$results"
export IPF_BATCH_KEY=sample_id IPF_CELLTYPE_KEY=cell_type IPF_THREADS="${SLURM_CPUS_PER_TASK:-32}"
export IPF_SEED=0 IPF_EPOCHS=500 IPF_BATCH_SIZE=32
export IPF_SUPERVISED_TARGET_KEY=cell_type IPF_SUPERVISED_TARGET_WEIGHTS="0.2 0.5 0.7 0.9"
export IPF_SUPERVISED_NEIGHBORS=15 IPF_SUPERVISED_MIN_DIST=0.5 MPLBACKEND=Agg
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python_exe=/home/lijia/luozhixiong/miniconda3/envs/methylvi/bin/python
cd "$project"
"$python_exe" Scripts/Methylvi/shared/04_train_methylvi.py
"$python_exe" Scripts/Methylvi/vmr/03_plot_vmr_methylvi_umap.py \
  --input "$results/methylvi_embedding.h5ad" --output-dir "$results"
"$python_exe" Scripts/Methylvi/shared/05_plot_supervised_umap.py
"$python_exe" Scripts/Methylvi/shared/07_plot_methylation_qc.py --input "$input" --output-dir "$results/qc"
touch "$route/model.COMPLETE"
