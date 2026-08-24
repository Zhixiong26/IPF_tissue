#!/usr/bin/env bash
# Run after the matching BAM array task; each sample has exactly two FASTQ partitions.
#SBATCH --job-name=ipf_qc_fastq
#SBATCH --partition=fat
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --output=/home/lijia/luozhixiong/IPF_tissue/Scripts/QC/logs/qc_fastq_%A_%a.out
#SBATCH --error=/home/lijia/luozhixiong/IPF_tissue/Scripts/QC/logs/qc_fastq_%A_%a.err

set -euo pipefail

project_dir=/home/lijia/luozhixiong/IPF_tissue
python_exe=/home/lijia/jiangyuanpei/miniforge3/envs/allcools/bin/python

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
whitelist=${sample_output}/02_bam.tsv
output_table=${sample_output}/01_fastq.tsv

if [[ ! -r "$whitelist" ]]; then
  printf 'Missing BAM-derived whitelist: %s\n' "$whitelist" >&2
  exit 3
fi
if [[ -e "$output_table" ]]; then
  printf 'Refusing to overwrite: %s\n' "$output_table" >&2
  exit 4
fi

fastq_args=()
while IFS= read -r r1; do
  r2=${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}
  if [[ ! -r "$r2" ]]; then
    printf 'Missing paired R2: %s\n' "$r2" >&2
    exit 5
  fi
  fastq_args+=(--r1 "$r1" --r2 "$r2")
done < <(find -L "${project_dir}/Data/Raw_fastq/Met/${batch}" -type f \
          -name '*_R1_001.fastq.gz' -print | sort)

if (( ${#fastq_args[@]} != 8 )); then
  printf 'Expected exactly two R1/R2 partitions for %s; argument count=%s\n' \
    "$batch" "${#fastq_args[@]}" >&2
  exit 6
fi

"$python_exe" "${project_dir}/Scripts/QC/01_FASTQ_per_cell_reads.py" \
  --sample-id "$sample_id" \
  "${fastq_args[@]}" \
  --whitelist "$whitelist" \
  --workers "${SLURM_CPUS_PER_TASK}" \
  --progress-interval 300 \
  --output "$output_table"

printf 'Completed FASTQ counting for %s: %s\n' "$sample_id" "$output_table"
