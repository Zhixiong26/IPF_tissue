#!/usr/bin/env bash
# One-click submit: common stages -> three parallel threshold branches -> summary.
# Usage: submit_methscan_pipeline.sh OUTPUT_DIR
set -Eeuo pipefail
script_dir=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methscan
if (( $# != 1 )); then echo "Usage: $0 OUTPUT_DIR" >&2; exit 2; fi
output_dir=$1
common_job=$(sbatch --parsable "$script_dir/run_methscan_common.sbatch" "$output_dir")
prepare_job=$(sbatch --parsable --dependency="afterok:${common_job}" "$script_dir/run_methscan_qc_stage.sbatch" "$output_dir" prepare)
filter_job=$(sbatch --parsable --dependency="afterok:${prepare_job}" "$script_dir/run_methscan_qc_stage.sbatch" "$output_dir" filter)
smooth_job=$(sbatch --parsable --dependency="afterok:${filter_job}" "$script_dir/run_methscan_qc_stage.sbatch" "$output_dir" smooth)
branch_jobs=()
for threshold in 0.01 0.02 0.05; do
  branch_jobs+=("$(sbatch --parsable --dependency="afterok:${smooth_job}" "$script_dir/run_methscan_branch.sbatch" "$output_dir" "$threshold")")
done
dependency=$(IFS=:; echo "afterok:${branch_jobs[*]}")
summary_job=$(sbatch --parsable --dependency="$dependency" "$script_dir/run_methscan_summary.sbatch" "$output_dir")
printf 'common=%s\nprepare=%s\nfilter=%s\nsmooth=%s\nbranches=%s\nsummary=%s\n' "$common_job" "$prepare_job" "$filter_job" "$smooth_job" "${branch_jobs[*]}" "$summary_job"
