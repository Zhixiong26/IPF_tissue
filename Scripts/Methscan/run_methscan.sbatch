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

set -euo pipefail

script_dir=/home/lijia/luozhixiong/IPF_tissue/Scripts/Methscan
# shellcheck source=00_methscan_config.sh
source "${script_dir}/00_methscan_config.sh"

if (( $# < 1 || $# > 3 )); then
  printf 'Usage: sbatch %s OUTPUT_DIR [QC_TABLE] [QC_COLUMN]\n' "$0" >&2
  exit 2
fi
output_dir=$1
qc_table=${2:-}
qc_column=${3:-pass_final_qc}

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
if [[ -n "$qc_table" && ! -r "$qc_table" ]]; then
  printf 'QC table is not readable: %s\n' "$qc_table" >&2
  exit 7
fi

available_kb=$(df -Pk "${project_dir}/Results/Methscan" | awk 'NR==2 {print $4}')
required_kb=$((min_free_gb * 1024 * 1024))
if (( available_kb < required_kb )); then
  printf 'At least %s GiB free is required; only %s KiB is available.\n' "$min_free_gb" "$available_kb" >&2
  exit 8
fi
mkdir -p "$output_dir"
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
  --max-cells "$max_cells"
  --require-index
)
if [[ -n "$qc_table" ]]; then
  manifest_args+=(--qc-table "$qc_table" --qc-column "$qc_column")
fi
"$methscan_python" "$script_dir/01_prepare_allc_inputs.py" "${manifest_args[@]}"

"$methscan_python" "$script_dir/02_prepare_methscan.py" \
  --methscan "$methscan_exe" \
  --manifest "$output_dir/00_manifest/input_manifest.tsv" \
  --chunksize "$prepare_chunksize" \
  --output-dir "$output_dir/01_prepared"

"$methscan_exe" filter \
  --min-sites "$min_sites" \
  --min-meth "$min_meth" \
  --max-meth "$max_meth" \
  "$output_dir/01_prepared" \
  "$output_dir/02_filtered"

"$methscan_exe" smooth \
  --bandwidth "$smooth_bandwidth" \
  "$output_dir/02_filtered"

mkdir -p "$output_dir/03_scan"
"$methscan_exe" scan \
  --bandwidth "$scan_bandwidth" \
  --stepsize "$scan_stepsize" \
  --var-threshold "$scan_var_threshold" \
  --min-cells "$scan_min_cells" \
  --threads "${SLURM_CPUS_PER_TASK:-$threads}" \
  "$output_dir/02_filtered" \
  "$output_dir/03_scan/VMRs.bed"

"$methscan_exe" matrix \
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
if [[ -n "$qc_table" ]]; then
  scanpy_args+=(--cell-metadata "$qc_table")
fi
"$scanpy_python" "$script_dir/03_vmr_scanpy.py" "${scanpy_args[@]}"

if [[ -n "$tss_bed" ]]; then
  if [[ ! -r "$tss_bed" ]]; then
    printf 'Configured TSS BED is unreadable: %s\n' "$tss_bed" >&2
    exit 9
  fi
  mkdir -p "$output_dir/06_profile"
  "$methscan_exe" profile \
    --strand-column "$tss_strand_column" \
    "$tss_bed" \
    "$output_dir/02_filtered" \
    "$output_dir/06_profile/TSS_profile.csv"
fi

"$scanpy_python" "$script_dir/04_summarize_run.py" --run-dir "$output_dir"

printf 'MethSCAn preliminary workflow completed: %s\n' "$output_dir"
