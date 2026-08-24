#!/usr/bin/env bash
# CYL/ZCP transcriptome E Scanpy workflow.
#SBATCH --job-name=ipf_scanpy_e
#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=/home/lijia/luozhixiong/IPF_tissue/Scripts/Scanpy/logs/scanpy_e_%j.out
#SBATCH --error=/home/lijia/luozhixiong/IPF_tissue/Scripts/Scanpy/logs/scanpy_e_%j.err

set -euo pipefail

project_dir=/home/lijia/luozhixiong/IPF_tissue
python_exe=/home/lijia/jiangyuanpei/miniforge3/envs/allcools/bin/python

if (( $# != 1 )); then
  printf 'Usage: sbatch %s OUTPUT_DIR\n' "$0" >&2
  exit 2
fi

output_dir=$1
if [[ "$output_dir" != "${project_dir}/Data/"* ]]; then
  printf 'OUTPUT_DIR must be under %s/Data: %s\n' "$project_dir" "$output_dir" >&2
  exit 3
fi
if [[ -e "$output_dir/rna_e_scanpy.h5ad" ]]; then
  printf 'Refusing to overwrite: %s\n' "$output_dir/rna_e_scanpy.h5ad" >&2
  exit 4
fi

"$python_exe" "$project_dir/Scripts/Scanpy/01_run_e_scanpy.py" \
  --cyl-zip "$project_dir/Data/Matrix/25100718_CYL_E/filtered_feature_bc_matrix.zip" \
  --zcp-zip "$project_dir/Data/Matrix/25100718_ZCP_E/filtered_feature_bc_matrix.zip" \
  --output-dir "$output_dir"
