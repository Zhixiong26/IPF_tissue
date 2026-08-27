#!/usr/bin/env bash
# Run selection and ALLC-to-COV only. Later QC stages are separate jobs.
# Usage: sbatch [resource options] run_methscan_common.sbatch OUTPUT_DIR
#SBATCH --job-name=ipf_ms_common
#SBATCH --partition=cpu
#SBATCH --nodelist=cu03
#SBATCH --nodes=1
#SBATCH --cpus-per-task=55
#SBATCH --mem=250G
#SBATCH --time=5-00:00:00
#SBATCH --output=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methscan/logs/methscan_common_%j.out
#SBATCH --error=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methscan/logs/methscan_common_%j.err

set -Eeuo pipefail
export PYTHONUNBUFFERED=1
script_dir=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methscan
source "$script_dir/00_methscan_config.sh"
if (( $# != 1 )); then echo "Usage: sbatch $0 OUTPUT_DIR" >&2; exit 2; fi
output_dir=$1
case "$output_dir" in "$project_dir/Results/Methscan/"*) ;; *) echo "OUTPUT_DIR must be under $project_dir/Results/Methscan/" >&2; exit 3;; esac
[[ ! -e "$output_dir" ]] || { echo "Output already exists: $output_dir" >&2; exit 4; }
[[ -x "$methscan_exe" && -x "$methscan_python" && -x "$scanpy_python" ]] || { echo "Configured executable missing" >&2; exit 5; }
[[ -d "$allc_source" && -r "$rna_annotation_table" ]] || { echo "Input path unavailable" >&2; exit 6; }
mkdir -p "${project_dir}/Results/Methscan"
available_kb=$(df -Pk "${project_dir}/Results/Methscan" | awk 'NR==2 {print $4}')
required_kb=$((min_free_gb * 1024 * 1024))
(( available_kb >= required_kb )) || { echo "Insufficient free space" >&2; exit 7; }
mkdir -p "$output_dir"
status_file="$output_dir/common_status.tsv"
printf 'stage\tstatus\ttimestamp\n' > "$status_file"
printf 'key\tvalue\nallc_source\t%s\nrna_annotation_table\t%s\nscan_var_thresholds\t%s\n' "$allc_source" "$rna_annotation_table" "${scan_var_thresholds[*]}" > "$output_dir/effective_configuration.tsv"
record_stage() { printf '%s\t%s\t%s\n' "$1" "$2" "$(date -Is)" >> "$status_file"; }
current_stage=initialization
on_error() { local code=$?; trap - ERR; record_stage "$current_stage" "failed(exit=$code)"; exit "$code"; }
trap on_error ERR
run_stage() { current_stage=$1; shift; record_stage "$current_stage" running; "$methscan_python" "$script_dir/06_run_with_resources.py" --output "$output_dir/${current_stage}.resources.json" -- "$@"; record_stage "$current_stage" complete; }
{
  "$methscan_exe" --version; "$methscan_python" --version
  "$scanpy_python" -c 'import scanpy; print("scanpy=" + scanpy.__version__)'
} > "$output_dir/software_versions.txt"
run_stage 00_scanpy_selection "$methscan_python" "$script_dir/01_select_scanpy_cells.py" --allc-source "$allc_source" --scanpy-annotation "$rna_annotation_table" --sample CYL --sample ZCP --exclude-cell-type "$rna_exclude_cell_type" --output-dir "$output_dir/00_scanpy_selected"
run_stage 01_convert_cov "$methscan_python" "$script_dir/02_convert_allc_to_cov.py" --manifest "$output_dir/00_scanpy_selected/input_manifest.tsv" --output-dir "$output_dir/01_cov" --workers "${SLURM_CPUS_PER_TASK:-$cov_conversion_workers}" --compresslevel "$cov_compresslevel"
touch "$output_dir/common.COMPLETE"
record_stage common workflow_complete
printf 'Common MethSCAn stages completed: %s\n' "$output_dir"
