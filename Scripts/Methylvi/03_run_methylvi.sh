#!/usr/bin/env bash
# Unified .cov -> ALLCools -> MethylVI entry point.
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$here/00_methylvi_config.sh"
stage=${1:-all}

# ALLCools validates helper programs (tabix/bgzip/bedtools) through PATH.
export PATH="$IPF_ALLCOOLS_ENV/bin:$PATH"
numeric_threads="${IPF_NUMERIC_THREADS:-1}"
export OMP_NUM_THREADS="$numeric_threads"
export MKL_NUM_THREADS="$numeric_threads"
export OPENBLAS_NUM_THREADS="$numeric_threads"
export NUMEXPR_NUM_THREADS="$numeric_threads"
export MPLBACKEND="${MPLBACKEND:-Agg}"

allc_python="$IPF_ALLCOOLS_ENV/bin/python"
allcools="$IPF_ALLCOOLS_ENV/bin/allcools"
mvi_python="$IPF_METHYLVI_ENV/bin/python"

require_file() {
  [[ -s "$1" ]] || { echo "ERROR: required file is missing or empty: $1" >&2; exit 1; }
}

verify() {
  [[ -d "$IPF_COV_DIR" ]] || { echo "ERROR: coverage directory missing: $IPF_COV_DIR" >&2; exit 1; }
  require_file "$IPF_ANNOTATION"
  require_file "$IPF_CHROM_SIZES"
  require_file "$IPF_BLACKLIST"
  [[ -x "$allc_python" && -x "$allcools" ]] || { echo "ERROR: invalid ALLCools environment: $IPF_ALLCOOLS_ENV" >&2; exit 1; }
  [[ -x "$mvi_python" ]] || { echo "ERROR: invalid MethylVI environment: $IPF_METHYLVI_ENV" >&2; exit 1; }
  command -v bedtools >/dev/null || { echo "ERROR: bedtools is required for blacklist filtering" >&2; exit 1; }
  command -v intersectBed >/dev/null || { echo "ERROR: intersectBed is required for blacklist filtering" >&2; exit 1; }
  "$allc_python" -c 'import ALLCools,anndata,pandas,scanpy; print("ALLCools environment OK", ALLCools.__version__)'
  "$mvi_python" -c 'import anndata,mudata,scanpy,scvi,torch; from scvi.external import METHYLVI; print("MethylVI environment OK", scvi.__version__)'
  "$allc_python" "$here/01_prepare_allcools.py" --verify-only
}

prepare_allc() {
  "$allc_python" "$here/01_prepare_allcools.py"
}

generate_mcds() {
  require_file "$IPF_ALLC_TABLE"
  mkdir -p "$IPF_ALLCOOLS_ROOT"
  table_hash=$(sha256sum "$IPF_ALLC_TABLE" | awk '{print $1}')
  config="table_sha256=$table_hash bin_size=$IPF_BIN_SIZE context=$IPF_MC_CONTEXT hypo_cutoff=$IPF_HYPO_SCORE_CUTOFF"
  config_file="$IPF_ALLCOOLS_ROOT/mcds.config.txt"
  complete="$IPF_ALLCOOLS_ROOT/mcds.COMPLETE"
  if [[ -e "$complete" ]]; then
    require_file "$config_file"
    [[ $(<"$config_file") == "$config" ]] || {
      echo "ERROR: existing MCDS was built with different inputs/parameters; use a new IPF_MVI_ROOT" >&2
      exit 1
    }
    echo "Existing compatible MCDS detected: $IPF_MCDS"
    return
  fi
  if [[ -e "$IPF_MCDS" ]]; then
    echo "ERROR: incomplete MCDS path exists without completion marker: $IPF_MCDS" >&2
    exit 1
  fi
  "$allcools" generate-dataset \
    --allc_table "$IPF_ALLC_TABLE" \
    --output_path "$IPF_MCDS" \
    --chrom_size_path "$IPF_CHROM_SIZES" \
    --obs_dim cell \
    --cpu "$IPF_THREADS" \
    --chunk_size 10 \
    --regions chrom5k "$IPF_BIN_SIZE" \
    --quantifiers chrom5k hypo-score "$IPF_MC_CONTEXT" "cutoff=$IPF_HYPO_SCORE_CUTOFF"
  printf '%s\n' "$config" > "$config_file"
  touch "$complete"
}

cluster_allcools() {
  require_file "$IPF_ALLCOOLS_ROOT/mcds.config.txt"
  annotation_hash=$(sha256sum "$IPF_ANNOTATION" | awk '{print $1}')
  mcds_config_hash=$(sha256sum "$IPF_ALLCOOLS_ROOT/mcds.config.txt" | awk '{print $1}')
  blacklist_hash=$(md5sum "$IPF_BLACKLIST" | awk '{print $1}')
  cluster_config="mcds_config_sha256=$mcds_config_hash annotation_sha256=$annotation_hash blacklist_md5=$blacklist_hash blacklist_fraction=$IPF_BLACKLIST_FRACTION bin_cutoff=$IPF_BINARIZE_CUTOFF hypo_percent=$IPF_HYPO_PERCENT lsi_components=$IPF_LSI_COMPONENTS p_cutoff=$IPF_LSI_P_CUTOFF neighbors=$IPF_ALLCOOLS_NEIGHBORS leiden=$IPF_ALLCOOLS_LEIDEN_RESOLUTION repeats=$IPF_CONSENSUS_LEIDEN_REPEATS consensus_leiden=$IPF_CONSENSUS_LEIDEN_RESOLUTION seed=$IPF_SEED"
  cluster_config_file="$IPF_ALLCOOLS_ROOT/cluster.config.txt"
  if [[ -s "$IPF_ALLCOOLS_H5AD" ]]; then
    require_file "$cluster_config_file"
    [[ $(<"$cluster_config_file") == "$cluster_config" ]] || {
      echo "ERROR: existing clustered H5AD uses different parameters; use a new IPF_MVI_ROOT" >&2
      exit 1
    }
    echo "Existing compatible ALLCools H5AD detected: $IPF_ALLCOOLS_H5AD"
    return
  fi
  "$allc_python" "$here/02_cluster_allcools.py"
  printf '%s\n' "$cluster_config" > "$cluster_config_file"
}

run_allcools() {
  prepare_allc
  generate_mcds
  cluster_allcools
}

build_methylvi() {
  require_file "$IPF_ALLCOOLS_H5AD"
  "$mvi_python" "$here/03_build_methylvi_from_allcools.py"
}

train_methylvi() {
  require_file "$IPF_MVI_INPUT"
  "$mvi_python" "$here/04_train_methylvi.py"
}

supervised_umap() {
  require_file "$IPF_MVI_RESULTS/methylvi_embedding.h5ad"
  "$mvi_python" "$here/05_plot_supervised_umap.py"
}

plots_before_methylvi() {
  require_file "$IPF_ALLCOOLS_H5AD"
  "$allc_python" "$here/06_plot_allcools_umap.py"
}

plots_after_methylvi() {
  require_file "$IPF_MVI_RESULTS/methylvi_embedding.h5ad"
  "$mvi_python" "$here/07_plot_methylvi_umap.py"
}

configure_smoke() {
  export IPF_MVI_ROOT="$IPF_PROJECT_DIR/Results/MethylVI_30wcov_allcools/smoke"
  export IPF_ALLCOOLS_ROOT="$IPF_MVI_ROOT/allcools_5kb"
  export IPF_ALLC_DIR="$IPF_ALLCOOLS_ROOT/input_allc"
  export IPF_ALLC_TABLE="$IPF_ALLCOOLS_ROOT/selected_cells.allc.tsv"
  export IPF_MCDS="$IPF_ALLCOOLS_ROOT/mcg_5kb.mcds"
  export IPF_ALLCOOLS_H5AD="$IPF_ALLCOOLS_ROOT/mcg_5kb.clustered.h5ad"
  export IPF_MVI_INPUT="$IPF_MVI_ROOT/ipf_allcools_5kb_methylvi_input.h5mu"
  export IPF_MVI_RESULTS="$IPF_MVI_ROOT/results"
  export IPF_MAX_CELLS="${IPF_SMOKE_CELLS:-100}"
  export IPF_BALANCED_COHORTS=1
  export IPF_EPOCHS="${IPF_SMOKE_EPOCHS:-2}"
}

case "$stage" in
  verify) verify ;;
  prepare) verify; prepare_allc; generate_mcds ;;
  cluster) cluster_allcools ;;
  allcools) verify; run_allcools ;;
  build) build_methylvi ;;
  train) train_methylvi ;;
  plots-before) plots_before_methylvi ;;
  plots-after) plots_after_methylvi ;;
  supervised) supervised_umap ;;
  all) verify; run_allcools; build_methylvi; train_methylvi; supervised_umap ;;
  smoke) configure_smoke; verify; run_allcools; build_methylvi; train_methylvi ;;
  *) echo "Usage: bash $0 {verify|prepare|cluster|allcools|build|train|plots-before|plots-after|supervised|smoke|all}" >&2; exit 2 ;;
esac
