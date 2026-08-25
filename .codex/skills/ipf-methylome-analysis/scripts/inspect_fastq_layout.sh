#!/usr/bin/env bash
# Read-only inventory of FASTQ-like files in the project-linked raw-data tree.
set -euo pipefail

project_dir="${1:-/home/lijia/luozhixiong/IPF_tissue}"
raw_dir="${project_dir}/Data/Raw_fastq"

if [[ ! -d "${raw_dir}" ]]; then
  printf 'Raw FASTQ directory not found: %s\n' "${raw_dir}" >&2
  exit 1
fi

printf 'group\tbatch\tsubset\tfilename\tbytes\n'
while IFS= read -r -d '' fastq; do
  filename="${fastq##*/}"
  subset_path="${fastq%/*}"
  subset="${subset_path##*/}"
  batch_path="${subset_path%/*}"
  batch="${batch_path##*/}"
  group_path="${batch_path%/*}"
  group="${group_path##*/}"
  bytes=$(stat -c '%s' "${fastq}")
  printf '%s\t%s\t%s\t%s\t%s\n' "${group}" "${batch}" "${subset}" "${filename}" "${bytes}"
done < <(find -L "${raw_dir}" -mindepth 4 -maxdepth 4 -type f \( -iname '*.fastq.gz' -o -iname '*.fq.gz' -o -iname '*.fastq' -o -iname '*.fq' \) -print0 | sort -z)
