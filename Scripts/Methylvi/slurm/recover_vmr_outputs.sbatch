#!/usr/bin/env bash
# Complete VMR models whose training finished before the plotting-wrapper fix.
# This is intentionally a post-training recovery stage; it never retrains.
#SBATCH --job-name=mvi_vmr_recover
#SBATCH --account=lijia
#SBATCH --partition=fat
#SBATCH --nodelist=fat01
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=45
#SBATCH --mem=90G
#SBATCH --time=2-00:00:00
#SBATCH --output=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/mvi_vmr_recover_%j.out
#SBATCH --error=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/mvi_vmr_recover_%j.err

set -Eeuo pipefail
project=/home/lijia/luozhixiong/IPF_tissue
root="$project/Results/MethylVI_clean6264_10k30k"
python_exe=/home/lijia/luozhixiong/miniconda3/envs/methylvi/bin/python
export IPF_BATCH_KEY=sample_id IPF_CELLTYPE_KEY=cell_type IPF_SEED=0
export IPF_THREADS=45 IPF_SUPERVISED_TARGET_KEY=cell_type
export IPF_SUPERVISED_TARGET_WEIGHTS="0.2 0.5 0.7 0.9"
export IPF_SUPERVISED_NEIGHBORS=15 IPF_SUPERVISED_MIN_DIST=0.5 MPLBACKEND=Agg
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

for version in 0.01 0.02 0.05; do
  for features in 10000 30000; do
    route="$root/vmr/var_${version}/features_${features}"
    input="$route/input.h5mu"
    results="$route/results"
    [[ -s "$input" && -s "$results/methylvi_embedding.h5ad" && -s "$results/run_summary.json" ]] || {
      echo "Missing completed training output: $route" >&2
      exit 1
    }
    export IPF_MVI_INPUT="$input" IPF_MVI_RESULTS="$results" VMR_MVI_RESULTS="$results"
    "$python_exe" Scripts/Methylvi/vmr/03_plot_vmr_methylvi_umap.py \
      --input "$results/methylvi_embedding.h5ad" --output-dir "$results"
    "$python_exe" Scripts/Methylvi/shared/05_plot_supervised_umap.py
    "$python_exe" Scripts/Methylvi/shared/07_plot_methylation_qc.py \
      --input "$input" --output-dir "$results/qc"
    touch "$route/model.COMPLETE"
  done
done
printf 'Recovered all six VMR post-training outputs under %s\n' "$root"
