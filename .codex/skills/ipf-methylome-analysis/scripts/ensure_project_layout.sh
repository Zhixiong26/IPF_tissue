#!/usr/bin/env bash
# Ensure the stable top-level directory contract used by this project skill.
set -euo pipefail

project_dir="${1:-/home/lijia/luozhixiong/IPF_tissue}"
stages=(Environment Scanpy Methscan Methylvi)
skill_script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

if (( $# > 1 )); then
  stages+=("${@:2}")
fi

if [[ ! -d "${project_dir}" ]]; then
  printf 'Project directory not found: %s\n' "${project_dir}" >&2
  exit 1
fi

for directory in Scripts Results Supplementary; do
  mkdir -p "${project_dir}/${directory}"
  printf '%s\n' "${project_dir}/${directory}"
done

for parent in Scripts Results; do
  for stage in "${stages[@]}"; do
    if [[ ! "${stage}" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
      printf 'Invalid stage directory name: %s\n' "${stage}" >&2
      exit 2
    fi
    mkdir -p "${project_dir}/${parent}/${stage}"
    printf '%s\n' "${project_dir}/${parent}/${stage}"
    if [[ "${parent}" == "Scripts" ]]; then
      mkdir -p "${project_dir}/${parent}/${stage}/logs"
      printf '%s\n' "${project_dir}/${parent}/${stage}/logs"
    fi
  done
done

python3 "${skill_script_dir}/init_stage_docs.py" "${project_dir}" "${stages[@]}"
