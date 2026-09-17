#!/usr/bin/env bash
# Submit pooled CYL+ZCP pairwise cell-type DMR analysis only.
# Usage: submit_methscan_pooled_dmr.sh OUTPUT_DIR
set -Eeuo pipefail
script_dir=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methscan
if (( $# != 1 )); then echo "Usage: $0 OUTPUT_DIR" >&2; exit 2; fi
output_dir=$1
diff_job=$(sbatch --parsable "$script_dir/run_methscan_pooled_methdiff.sbatch" "$output_dir")
printf 'pooled_methdiff=%s\n' "$diff_job"
