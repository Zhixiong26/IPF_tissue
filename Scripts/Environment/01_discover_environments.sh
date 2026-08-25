#!/usr/bin/env bash
# Print a read-only TSV inventory of active tools and known Conda environments.
set -euo pipefail

printf 'record_type\tname\tpath\tversion\tavailable_tools\n'

tools=(python3 conda mamba samtools bismark bismark_methylation_extractor bowtie2 fastqc multiqc methscan allcools bgzip tabix)
for tool in "${tools[@]}"; do
  tool_path=$(command -v "${tool}" 2>/dev/null || true)
  if [[ -z "${tool_path}" ]]; then
    printf 'active_tool\t%s\t\tmissing\t\n' "${tool}"
    continue
  fi
  tool_version=$("${tool_path}" --version 2>&1 | sed -n '1p' || true)
  tool_version=${tool_version//$'\t'/ }
  printf 'active_tool\t%s\t%s\t%s\t\n' "${tool}" "${tool_path}" "${tool_version:-unknown}"
done

env_roots=(
  /home/lijia/luozhixiong/miniconda3/envs
  /home/lijia/jiangyuanpei/miniforge3/envs
)

for env_root in "${env_roots[@]}"; do
  [[ -d "${env_root}" ]] || continue
  for env_path in "${env_root}"/*; do
    [[ -d "${env_path}" ]] || continue
    env_name=${env_path##*/}
    python_version=none
    if [[ -x "${env_path}/bin/python" ]]; then
      python_version=$("${env_path}/bin/python" --version 2>&1 | sed -n '1p' || true)
    fi
    available=""
    for tool in "${tools[@]}"; do
      if [[ -x "${env_path}/bin/${tool}" ]]; then
        available="${available}${available:+,}${tool}"
      fi
    done
    printf 'conda_env\t%s\t%s\t%s\t%s\n' "${env_name}" "${env_path}" "${python_version}" "${available}"
  done
done
