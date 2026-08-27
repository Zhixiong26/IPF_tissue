#!/usr/bin/env bash
# Independent VMR-region MethylVI configuration.

export VMR_PROJECT_DIR="${VMR_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
export VMR_SCRIPT_DIR="${VMR_SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
export VMR_COV_DIR="${VMR_COV_DIR:-${VMR_PROJECT_DIR}/Data/30wcov}"
export VMR_EXISTING_ALLC_DIR="${VMR_EXISTING_ALLC_DIR:-/home/lijia/jiangyuanpei/methscan/xunyin/IPF_tissue/allcools_5kbin/input_allc}"
# The current formal MethSCAn run supplies the selected, original ALLC paths.
# Its VMR BED becomes available only after the upstream scan stage completes.
export VMR_METHSCAN_RUN_DIR="${VMR_METHSCAN_RUN_DIR:-${VMR_PROJECT_DIR}/Results/Methscan/CYL_ZCP_full_20260826_final}"
export VMR_METHSCAN_VARIANCE="${VMR_METHSCAN_VARIANCE:-0.01}"
export VMR_INPUT_MANIFEST="${VMR_INPUT_MANIFEST:-${VMR_METHSCAN_RUN_DIR}/00_scanpy_selected/input_manifest.tsv}"
# Set this to the selected MethSCAn branch, for example:
# Results/Methscan/<run>/04_scan/var_0.01/VMRs.bed
export VMR_SOURCE_BED="${VMR_SOURCE_BED:-${IPF_METHSCAN_VMR_SOURCE:-${VMR_METHSCAN_RUN_DIR}/04_scan/var_${VMR_METHSCAN_VARIANCE}/VMRs.bed}}"
export VMR_CHROM_SIZES="${VMR_CHROM_SIZES:-${VMR_PROJECT_DIR}/Supplementary/hg38.canonical.chrom.sizes}"
export VMR_BLACKLIST="${VMR_BLACKLIST:-${VMR_PROJECT_DIR}/Supplementary/ENCFF356LFX_GRCh38_blacklist.bed.gz}"
export VMR_BLACKLIST_MD5="${VMR_BLACKLIST_MD5:-393688b4f06c9ce26165d47433dd8c37}"
export VMR_BLACKLIST_FRACTION="${VMR_BLACKLIST_FRACTION:-0.2}"
export VMR_ANNOTATION="${VMR_ANNOTATION:-${VMR_PROJECT_DIR}/Supplementary/manual_celltype_annotation.tsv}"

export VMR_RESULTS_ROOT="${VMR_RESULTS_ROOT:-${VMR_PROJECT_DIR}/Results/MethylVI_30wcov_vmrs_blacklist_f0p2}"
export VMR_INPUT_DIR="${VMR_INPUT_DIR:-${VMR_RESULTS_ROOT}/input}"
export VMR_FILTERED_BED="${VMR_FILTERED_BED:-${VMR_INPUT_DIR}/vmr_canonical_blacklist_f0p2.bed}"
export VMR_ALLC_TABLE="${VMR_ALLC_TABLE:-${VMR_INPUT_DIR}/selected_cells.allc.tsv}"
export VMR_MVI_INPUT="${VMR_MVI_INPUT:-${VMR_RESULTS_ROOT}/ipf_vmr_methylvi_input.h5mu}"
export VMR_COUNT_ROWS="${VMR_COUNT_ROWS:-${VMR_RESULTS_ROOT}/count_rows}"
export VMR_MVI_RESULTS="${VMR_MVI_RESULTS:-${VMR_RESULTS_ROOT}/results}"

export VMR_METHYLVI_ENV="${VMR_METHYLVI_ENV:-/home/lijia/luozhixiong/miniconda3/envs/methylvi}"
# Set a positive value only for the legacy coverage-directory fallback.
export VMR_EXPECTED_CELLS="${VMR_EXPECTED_CELLS:-0}"
# Retain a VMR when coverage is observed in more than 200 cells (~3.06%).
export VMR_MIN_COVERED_CELLS="${VMR_MIN_COVERED_CELLS:-200}"
export VMR_THREADS="${VMR_THREADS:-32}"
export VMR_SEED="${VMR_SEED:-0}"
export VMR_EPOCHS="${VMR_EPOCHS:-500}"
export VMR_BATCH_SIZE="${VMR_BATCH_SIZE:-32}"
export VMR_BATCH_KEY="${VMR_BATCH_KEY:-cohort}"
export VMR_CELLTYPE_KEY="${VMR_CELLTYPE_KEY:-manual_celltype}"

# Compatibility variables consumed by the shared MethylVI training/UMAP code.
export IPF_MVI_INPUT="$VMR_MVI_INPUT"
export IPF_MVI_RESULTS="$VMR_MVI_RESULTS"
export IPF_THREADS="$VMR_THREADS"
export IPF_SEED="$VMR_SEED"
export IPF_EPOCHS="$VMR_EPOCHS"
export IPF_BATCH_SIZE="$VMR_BATCH_SIZE"
export IPF_BATCH_KEY="$VMR_BATCH_KEY"
export IPF_CELLTYPE_KEY="$VMR_CELLTYPE_KEY"
export IPF_SUPERVISED_TARGET_KEY="${IPF_SUPERVISED_TARGET_KEY:-$VMR_CELLTYPE_KEY}"
export IPF_SUPERVISED_TARGET_WEIGHTS="${IPF_SUPERVISED_TARGET_WEIGHTS:-0.2 0.5 0.7 0.9}"
export IPF_SUPERVISED_NEIGHBORS="${IPF_SUPERVISED_NEIGHBORS:-15}"
export IPF_SUPERVISED_MIN_DIST="${IPF_SUPERVISED_MIN_DIST:-0.5}"
