#!/usr/bin/env bash
# Submit one new output directory for the complete three-threshold workflow.
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
source "$script_dir/00_methscan_config.sh"

if (( $# != 1 )); then echo "Usage: sbatch $0 OUTPUT_DIR" >&2; exit 2; fi
output_dir=$1
case "$output_dir" in
  "$project_dir/Results/Methscan/"*) ;;
  *) echo "OUTPUT_DIR must be under $project_dir/Results/Methscan/" >&2; exit 3;;
esac
[[ ! -e "$output_dir" ]] || { echo "Output already exists: $output_dir" >&2; exit 4; }
[[ -x "$methscan_exe" && -x "$methscan_python" && -x "$scanpy_python" ]] || { echo "Configured executable missing" >&2; exit 5; }
[[ -d "$allc_source" && -r "$rna_annotation_table" ]] || { echo "Input path unavailable" >&2; exit 6; }

mkdir -p "${project_dir}/Results/Methscan"
available_kb=$(df -Pk "${project_dir}/Results/Methscan" | awk 'NR==2 {print $4}')
required_kb=$((min_free_gb * 1024 * 1024))
(( available_kb >= required_kb )) || { echo "Insufficient free space" >&2; exit 7; }
mkdir -p "$output_dir"
stage_status="$output_dir/stage_status.tsv"
printf 'stage\tstatus\ttimestamp\n' > "$stage_status"
printf 'key\tvalue\nallc_source\t%s\nrna_annotation_table\t%s\nscan_var_thresholds\t%s\n' \
  "$allc_source" "$rna_annotation_table" "${scan_var_thresholds[*]}" > "$output_dir/effective_configuration.tsv"

current_stage=initialization
record_stage() { printf '%s\t%s\t%s\n' "$1" "$2" "$(date -Is)" >> "$stage_status"; }
on_error() { local code=$?; trap - ERR; record_stage "$current_stage" "failed(exit=$code)"; exit "$code"; }
on_signal() { local sig=$1; record_stage "$current_stage" "cancelled(signal=$sig)"; trap - "$sig"; kill -s "$sig" "$$"; }
trap on_error ERR
trap 'on_signal TERM' TERM
trap 'on_signal INT' INT
run_stage() {
  current_stage=$1; shift; record_stage "$current_stage" running
  "$methscan_python" "$script_dir/06_run_with_resources.py" --output "$output_dir/${current_stage}.resources.json" -- "$@"
  record_stage "$current_stage" complete
}

{
  "$methscan_exe" --version; "$methscan_python" --version
  "$scanpy_python" -c 'import scanpy; print("scanpy=" + scanpy.__version__)'
} > "$output_dir/software_versions.txt"

run_stage 00_scanpy_selection "$methscan_python" "$script_dir/01_select_scanpy_cells.py" \
  --allc-source "$allc_source" --scanpy-annotation "$rna_annotation_table" \
  --sample CYL --sample ZCP --exclude-cell-type "$rna_exclude_cell_type" \
  --output-dir "$output_dir/00_scanpy_selected"
run_stage 01_convert_cov "$methscan_python" "$script_dir/02_convert_allc_to_cov.py" \
  --manifest "$output_dir/00_scanpy_selected/input_manifest.tsv" \
  --output-dir "$output_dir/01_cov" \
  --workers "${SLURM_CPUS_PER_TASK:-$cov_conversion_workers}" \
  --compresslevel "$cov_compresslevel"
run_stage 02_prepare "$methscan_python" "$script_dir/03_prepare_methscan.py" \
  --methscan "$methscan_exe" --manifest "$output_dir/01_cov/input_manifest.tsv" \
  --chunksize "$prepare_chunksize" --output-dir "$output_dir/02_prepared"
run_stage 03_filter "$methscan_exe" filter --min-sites "$min_sites" --min-meth "$min_meth" \
  --max-meth "$max_meth" "$output_dir/02_prepared" "$output_dir/03_filtered"
run_stage 04_smooth "$methscan_exe" smooth --bandwidth "$smooth_bandwidth" "$output_dir/03_filtered"

mkdir -p "$output_dir/04_scan" "$output_dir/05_matrix" "$output_dir/06_scanpy"
for var_threshold in "${scan_var_thresholds[@]}"; do
  label="var_${var_threshold}"
  scan_dir="$output_dir/04_scan/$label"
  matrix_dir="$output_dir/05_matrix/$label"
  scanpy_dir="$output_dir/06_scanpy/$label"
  mkdir -p "$scan_dir"
  run_stage "05_scan_${label}" "$methscan_exe" scan --bandwidth "$scan_bandwidth" \
    --stepsize "$scan_stepsize" --var-threshold "$var_threshold" --min-cells "$scan_min_cells" \
    --threads "${SLURM_CPUS_PER_TASK:-$threads}" "$output_dir/03_filtered" "$scan_dir/VMRs.bed"
  run_stage "06_matrix_${label}" "$methscan_exe" matrix \
    --threads "${SLURM_CPUS_PER_TASK:-$threads}" "$scan_dir/VMRs.bed" \
    "$output_dir/03_filtered" "$matrix_dir"
  run_stage "07_scanpy_${label}" "$scanpy_python" "$script_dir/04_vmr_scanpy.py" \
    --matrix "$matrix_dir/mean_shrunken_residuals.csv.gz" \
    --expected-cell-ids "$output_dir/03_filtered/column_header.txt" \
    --cell-metadata "$output_dir/00_scanpy_selected/input_manifest.tsv" \
    --output-dir "$scanpy_dir" --min-region-cell-fraction "$vmr_min_cell_fraction" \
    --min-cell-regions "$cell_min_regions" --n-pcs "$pca_components" \
    --neighbors "$neighbors" --leiden-resolution "$leiden_resolution" --seed "$random_seed"
done
summary_args=(--run-dir "$output_dir")
for var_threshold in "${scan_var_thresholds[@]}"; do summary_args+=(--threshold "$var_threshold"); done
run_stage 08_summary "$scanpy_python" "$script_dir/05_summarize_run.py" "${summary_args[@]}"
current_stage=workflow; record_stage workflow complete
trap - ERR TERM INT
printf 'MethSCAn workflow completed: %s\n' "$output_dir"
