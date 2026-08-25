#!/usr/bin/env bash
# Preliminary CYL/ZCP MethSCAn VMR workflow; no job is submitted by this file alone.
#SBATCH --job-name=ipf_methscan
#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=5-00:00:00
#SBATCH --output=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methscan/logs/methscan_%j.out
#SBATCH --error=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methscan/logs/methscan_%j.err

set -Eeuo pipefail
export PYTHONUNBUFFERED=1

script_dir=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methscan
# shellcheck source=00_methscan_config.sh
source "${script_dir}/00_methscan_config.sh"

if (( $# != 1 )); then
  printf 'Usage: sbatch %s OUTPUT_DIR\n' "$0" >&2
  exit 2
fi
output_dir=$1

if [[ "$output_dir" != "${project_dir}/Results/Methscan/"* ]]; then
  printf 'OUTPUT_DIR must be a new path under %s/Results/Methscan/: %s\n' "$project_dir" "$output_dir" >&2
  exit 3
fi
if [[ -e "$output_dir" ]]; then
  printf 'Refusing to use an existing output path: %s\n' "$output_dir" >&2
  exit 4
fi
if [[ ! -x "$methscan_exe" || ! -x "$methscan_python" || ! -x "$scanpy_python" ]]; then
  printf 'A configured executable is missing or not executable.\n' >&2
  exit 5
fi
if [[ ! -e "$allc_source" ]]; then
  printf 'ALLC source is unavailable: %s\n' "$allc_source" >&2
  exit 6
fi
if [[ ! -r "$rna_annotation_table" ]]; then
  printf 'RNA annotation table is unreadable: %s\n' "$rna_annotation_table" >&2
  exit 7
fi

available_kb=$(df -Pk "${project_dir}/Results/Methscan" | awk 'NR==2 {print $4}')
required_kb=$((min_free_gb * 1024 * 1024))
if (( available_kb < required_kb )); then
  printf 'At least %s GiB free is required; only %s KiB is available.\n' "$min_free_gb" "$available_kb" >&2
  exit 8
fi
mkdir -p "$output_dir"
stage_status="$output_dir/stage_status.tsv"
printf 'stage\tstatus\ttimestamp\n' > "$stage_status"
current_stage=initialization

record_stage() {
  printf '%s\t%s\t%s\n' "$1" "$2" "$(date -Is)" >> "$stage_status"
}

on_error() {
  exit_code=$?
  trap - ERR
  record_stage "$current_stage" "failed(exit=${exit_code})"
  exit "$exit_code"
}
trap on_error ERR

on_signal() {
  signal_name=$1
  record_stage "$current_stage" "cancelled(signal=${signal_name})"
  trap - "$signal_name"
  kill -s "$signal_name" "$$"
}
trap 'on_signal TERM' TERM
trap 'on_signal INT' INT

run_stage() {
  current_stage=$1
  shift
  record_stage "$current_stage" running
  "$methscan_python" "$script_dir/05_run_with_resources.py" \
    --output "$output_dir/${current_stage}.resources.json" -- "$@"
  record_stage "$current_stage" complete
}

{
  "$methscan_exe" --version
  "$methscan_python" --version
  "$scanpy_python" -c 'import scanpy; print("scanpy=" + scanpy.__version__)'
} > "$output_dir/software_versions.txt"

manifest_args=(
  --allc-source "$allc_source"
  --output-dir "$output_dir/00_manifest"
  --sample CYL
  --sample ZCP
  --validation-records "$allc_validation_records"
  --workers "$allc_intake_workers"
  --max-cells 0
  --require-index
)
run_stage 00_manifest "$methscan_python" "$script_dir/01_prepare_allc_inputs.py" "${manifest_args[@]}"

run_stage 00_scanpy_selection "$methscan_python" "$script_dir/01_select_scanpy_cells.py" \
  --manifest "$output_dir/00_manifest/input_manifest.tsv" \
  --scanpy-annotation "$rna_annotation_table" \
  --exclude-cell-type "$rna_exclude_cell_type" \
  --max-cells "$max_cells" \
  --output-dir "$output_dir/00_scanpy_selected"

run_stage 01_prepare "$methscan_python" "$script_dir/02_prepare_methscan.py" \
  --methscan "$methscan_exe" \
  --manifest "$output_dir/00_scanpy_selected/input_manifest.tsv" \
  --chunksize "$prepare_chunksize" \
  --output-dir "$output_dir/01_prepared"

run_stage 02_filter "$methscan_exe" filter \
  --min-sites "$min_sites" \
  --min-meth "$min_meth" \
  --max-meth "$max_meth" \
  "$output_dir/01_prepared" \
  "$output_dir/02_filtered"

run_stage 03_smooth "$methscan_exe" smooth \
  --bandwidth "$smooth_bandwidth" \
  "$output_dir/02_filtered"

mkdir -p "$output_dir/03_scan"
run_stage 04_scan "$methscan_exe" scan \
  --bandwidth "$scan_bandwidth" \
  --stepsize "$scan_stepsize" \
  --var-threshold "$scan_var_threshold" \
  --min-cells "$scan_min_cells" \
  --threads "${SLURM_CPUS_PER_TASK:-$threads}" \
  "$output_dir/02_filtered" \
  "$output_dir/03_scan/VMRs.bed"

run_stage 05_matrix "$methscan_exe" matrix \
  --threads "${SLURM_CPUS_PER_TASK:-$threads}" \
  "$output_dir/03_scan/VMRs.bed" \
  "$output_dir/02_filtered" \
  "$output_dir/04_matrix"

scanpy_args=(
  --matrix "$output_dir/04_matrix/mean_shrunken_residuals.csv.gz" \
  --output-dir "$output_dir/05_scanpy" \
  --min-region-cell-fraction "$vmr_min_cell_fraction" \
  --min-cell-regions "$cell_min_regions" \
  --n-pcs "$pca_components" \
  --neighbors "$neighbors" \
  --leiden-resolution "$leiden_resolution" \
  --seed "$random_seed"
)
# RNA labels are fixed once at entry from the canonical Scanpy whitelist.
# MethSCAn filtering only removes cells; Scanpy inherits labels for retained IDs.
scanpy_args+=(--cell-metadata "$output_dir/00_scanpy_selected/input_manifest.tsv")
run_stage 06_scanpy "$scanpy_python" "$script_dir/03_vmr_scanpy.py" "${scanpy_args[@]}"

if [[ -n "$tss_bed" ]]; then
  if [[ ! -r "$tss_bed" ]]; then
    printf 'Configured TSS BED is unreadable: %s\n' "$tss_bed" >&2
    exit 9
  fi
  mkdir -p "$output_dir/06_profile"
  run_stage 07_profile "$methscan_exe" profile \
    --strand-column "$tss_strand_column" \
    "$tss_bed" \
    "$output_dir/02_filtered" \
    "$output_dir/06_profile/TSS_profile.csv"
fi

run_stage 08_summary "$scanpy_python" "$script_dir/04_summarize_run.py" --run-dir "$output_dir"

current_stage=workflow
record_stage workflow complete
trap - ERR TERM INT
printf 'MethSCAn preliminary workflow completed: %s\n' "$output_dir"
