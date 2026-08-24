#!/usr/bin/env bash
# CYL and ZCP run as an array; each sample scans its 156 prefix BAMs in parallel.
#SBATCH --job-name=ipf_qc_bam
#SBATCH --partition=fat
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=3-00:00:00
#SBATCH --output=/home/lijia/luozhixiong/IPF_tissue/Scripts/QC/logs/qc_bam_%A_%a.out
#SBATCH --error=/home/lijia/luozhixiong/IPF_tissue/Scripts/QC/logs/qc_bam_%A_%a.err

set -euo pipefail

project_dir=/home/lijia/luozhixiong/IPF_tissue
python_exe=/home/lijia/jiangyuanpei/miniforge3/envs/allcools/bin/python
samtools_exe=/home/lijia/jiangyuanpei/miniforge3/envs/allcools/bin/samtools

if (( $# != 1 )); then
  printf 'Usage: sbatch %s OUTPUT_ROOT\n' "$0" >&2
  exit 2
fi

output_root=$1
samples=(CYL ZCP)
batches=(25100718_CYL_Met 25100718_ZCP_Met)
array_index=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
sample_id=${samples[$array_index]}
batch=${batches[$array_index]}
sample_output=${output_root}/${sample_id}
output_table=${sample_output}/02_bam.tsv

if [[ -e "$output_table" ]]; then
  printf 'Refusing to overwrite: %s\n' "$output_table" >&2
  exit 3
fi
mkdir -p "$sample_output"

"$python_exe" "${project_dir}/Scripts/QC/02_BAM_per_cell_mapping.py" \
  --sample-id "$sample_id" \
  --bam-dir "${project_dir}/Data/Bam/${batch}/bismark" \
  --samtools "$samtools_exe" \
  --workers "${SLURM_CPUS_PER_TASK}" \
  --output "$output_table"

printf 'Completed BAM counting for %s: %s\n' "$sample_id" "$output_table"
