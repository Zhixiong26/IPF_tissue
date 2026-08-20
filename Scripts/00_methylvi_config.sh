#!/usr/bin/env bash
# IPF .cov -> MethylVI configuration.  Override any value before sourcing.

export IPF_PROJECT_DIR="${IPF_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export IPF_COV_DIR="${IPF_COV_DIR:-${IPF_PROJECT_DIR}/Data/30wcov}"
export IPF_ANNOTATION="${IPF_ANNOTATION:-${IPF_PROJECT_DIR}/Supplementary/manual_celltype_annotation.tsv}"
export IPF_CHROM_SIZES="${IPF_CHROM_SIZES:-${IPF_PROJECT_DIR}/Supplementary/hg38.canonical.chrom.sizes}"
export IPF_MVI_ROOT="${IPF_MVI_ROOT:-${IPF_PROJECT_DIR}/Results/MethylVI_30wcov}"
export IPF_MVI_INPUT="${IPF_MVI_INPUT:-${IPF_MVI_ROOT}/ipf_5kb_methylvi_input.h5mu}"
export IPF_MVI_RESULTS="${IPF_MVI_RESULTS:-${IPF_MVI_ROOT}/results}"
# ENCODE GRCh38 exclusion list from the reference MethylVI project.
export IPF_BLACKLIST="${IPF_BLACKLIST:-${IPF_PROJECT_DIR}/Supplementary/ENCFF356LFX_GRCh38_blacklist.bed.gz}"
export IPF_BLACKLIST_FRACTION="${IPF_BLACKLIST_FRACTION:-0.2}"

# The filename prefixes (CYL/ZCP) are the only available sample labels in the
# supplied directory.  They are used as the MethylVI batch covariate; replace
# this table/logic if donor IDs are available.
export IPF_BATCH_KEY="${IPF_BATCH_KEY:-cohort}"
export IPF_CELLTYPE_KEY="${IPF_CELLTYPE_KEY:-manual_celltype}"
export IPF_BIN_SIZE="${IPF_BIN_SIZE:-5000}"
export IPF_MIN_CELLS_PER_BIN="${IPF_MIN_CELLS_PER_BIN:-50}"
export IPF_MIN_COV_PER_BIN="${IPF_MIN_COV_PER_BIN:-20}"
export IPF_TOP_BINS="${IPF_TOP_BINS:-50000}"
export IPF_THREADS="${IPF_THREADS:-16}"
export IPF_SEED="${IPF_SEED:-0}"
export IPF_CONDA_ENV="${IPF_CONDA_ENV:-methylvi}"
