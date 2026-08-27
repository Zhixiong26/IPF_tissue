#!/usr/bin/env bash
# Validate and summarize all completed threshold branches.
# Usage: sbatch --dependency=afterok:JOB1:JOB2:JOB3 run_methscan_summary.sbatch OUTPUT_DIR
#SBATCH --job-name=ipf_ms_summary
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methscan/logs/methscan_summary_%j.out
#SBATCH --error=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methscan/logs/methscan_summary_%j.err

set -Eeuo pipefail
export PYTHONUNBUFFERED=1
script_dir=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methscan
source "$script_dir/00_methscan_config.sh"
if (( $# != 1 )); then echo "Usage: sbatch $0 OUTPUT_DIR" >&2; exit 2; fi
output_dir=$1
for threshold in "${scan_var_thresholds[@]}"; do [[ -f "$output_dir/var_${threshold}.COMPLETE" ]] || { echo "Missing completed branch: $threshold" >&2; exit 3; }; done
args=(--run-dir "$output_dir")
for threshold in "${scan_var_thresholds[@]}"; do args+=(--threshold "$threshold"); done
"$scanpy_python" "$script_dir/05_summarize_run.py" "${args[@]}"
touch "$output_dir/summary.COMPLETE"
printf 'MethSCAn summary completed: %s\n' "$output_dir"
