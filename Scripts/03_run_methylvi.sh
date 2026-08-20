#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/00_methylvi_config.sh"
stage=${1:-all}
if [[ ${IPF_SKIP_CONDA:-0} != 1 ]]; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$IPF_CONDA_ENV"
fi
python3 -c 'import sys; assert sys.version_info >= (3, 9), "MethylVI pipeline requires Python >= 3.9"'
case "$stage" in
  build) python3 "$here/01_build_cov_methylvi_input.py" ;;
  train) python3 "$here/02_train_methylvi.py" ;;
  all) python3 "$here/01_build_cov_methylvi_input.py" && python3 "$here/02_train_methylvi.py" ;;
  *) echo "Usage: bash $0 {build|train|all}" >&2; exit 2 ;;
esac
