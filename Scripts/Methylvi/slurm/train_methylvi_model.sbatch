#!/usr/bin/env bash
# Train and plot one of the eight method/version/feature models.
# Usage: sbatch train_methylvi_model.sbatch allcools none 10000
#        sbatch train_methylvi_model.sbatch vmr 0.01 30000
#SBATCH --job-name=mvi_train
#SBATCH --account=lijia
#SBATCH --partition=fat
#SBATCH --nodelist=fat01
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=45
#SBATCH --mem=90G
#SBATCH --time=7-00:00:00
#SBATCH --output=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/mvi_train_%j.out
#SBATCH --error=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/logs/mvi_train_%j.err

set -Eeuo pipefail
if (( $# != 3 )); then echo "Usage: sbatch $0 allcools|vmr VERSION 10000|30000" >&2; exit 2; fi
method=$1
version=$2
features=$3
case "$method" in allcools|vmr) ;; *) echo "Unsupported method: $method" >&2; exit 3;; esac
case "$features" in 10000|30000) ;; *) echo "Unsupported feature target: $features" >&2; exit 4;; esac
if [[ "$method" == vmr ]]; then
  case "$version" in 0.01|0.02|0.05) ;; *) echo "Unsupported VMR threshold: $version" >&2; exit 5;; esac
fi
project=/home/lijia/luozhixiong/IPF_tissue
source "$project/Scripts/Methylvi/00_experiment_config.sh"
if [[ "$method" == allcools ]]; then
  route="$MVI_EXPERIMENT_ROOT/allcools/features_${features}"
else
  route="$MVI_EXPERIMENT_ROOT/vmr/var_${version}/features_${features}"
fi
input="$route/input.h5mu"
results="$route/results"
[[ -s "$input" ]] || { echo "Missing input: $input" >&2; exit 6; }
export IPF_MVI_INPUT="$input"
export IPF_MVI_RESULTS="$results"
# The VMR plotting entry point reads its default output location from this
# route-specific variable; keep it aligned with the shared MethylVI results.
export VMR_MVI_RESULTS="$results"
export IPF_BATCH_KEY=sample_id
export IPF_CELLTYPE_KEY=cell_type
export IPF_THREADS="${SLURM_CPUS_PER_TASK:-45}"
export IPF_SEED=0
export IPF_EPOCHS=500
export IPF_BATCH_SIZE=32
export IPF_SUPERVISED_TARGET_KEY=cell_type
export IPF_SUPERVISED_TARGET_WEIGHTS="0.2 0.5 0.7 0.9"
export IPF_SUPERVISED_NEIGHBORS=15
export IPF_SUPERVISED_MIN_DIST=0.5
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd "$project"
python_exe="$MVI_METHYLVI_ENV/bin/python"
"$python_exe" - "$input" "$features" <<'PY'
import sys, mudata, numpy as np
x = mudata.read_h5mu(sys.argv[1], backed="r")["mCG"]
expected = int(sys.argv[2])
if x.n_vars != expected or x.n_obs != 6264:
    raise SystemExit(f"Input shape mismatch: {x.shape}, expected (6264, {expected})")
missing = {"sample_id", "condition", "cell_type"} - set(x.obs.columns)
if missing:
    raise SystemExit(f"Missing metadata: {sorted(missing)}")
if set(np.asarray(x.var["selection_rank"], dtype=int)) != set(range(1, expected + 1)):
    raise SystemExit("Invalid selection_rank")
print(f"Input preflight passed: {x.n_obs} cells x {x.n_vars} features")
PY
"$python_exe" Scripts/Methylvi/shared/04_train_methylvi.py
if [[ "$method" == allcools ]]; then
  "$python_exe" Scripts/Methylvi/allcools/07_plot_methylvi_umap.py --input "$results/methylvi_embedding.h5ad" --output-dir "$results"
else
  "$python_exe" Scripts/Methylvi/vmr/03_plot_vmr_methylvi_umap.py --input "$results/methylvi_embedding.h5ad" --output-dir "$results"
fi
"$python_exe" Scripts/Methylvi/shared/05_plot_supervised_umap.py
"$python_exe" Scripts/Methylvi/shared/07_plot_methylation_qc.py --input "$input" --output-dir "$results/qc"
touch "$route/model.COMPLETE"
