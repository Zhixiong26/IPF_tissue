#!/usr/bin/env bash
# Small compute-node check for the configured ALLC source and environments.
#SBATCH --job-name=ipf_methscan_check
#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methscan/logs/methscan_preflight_%j.out
#SBATCH --error=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methscan/logs/methscan_preflight_%j.err

set -euo pipefail

script_dir=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methscan
# shellcheck source=00_methscan_config.sh
source "${script_dir}/00_methscan_config.sh"

test -r "$allc_source"
test -x "$methscan_exe"
test -x "$methscan_python"
test -x "$scanpy_python"

preflight_dir=/tmp/ipf_methscan_preflight_${SLURM_JOB_ID:?SLURM_JOB_ID is required}
if [[ -e "$preflight_dir" ]]; then
  printf 'Temporary preflight path already exists: %s\n' "$preflight_dir" >&2
  exit 2
fi

"$methscan_python" "$script_dir/01_prepare_allc_inputs.py" \
  --allc-source "$allc_source" \
  --output-dir "$preflight_dir" \
  --sample CYL \
  --sample ZCP \
  --validation-records 100 \
  --max-cells 20 \
  --workers "${SLURM_CPUS_PER_TASK}" \
  --require-index \
  --verify-only

"$methscan_exe" --version
"$scanpy_python" -c 'import scanpy; print("scanpy=" + scanpy.__version__)'
printf 'MethSCAn compute-node preflight passed.\n'
