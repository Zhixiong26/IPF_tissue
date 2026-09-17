#!/usr/bin/env bash
# IPF .cov -> ALLCools mCG 5-kb clustering -> MethylVI configuration.
# Export a variable before sourcing this file to override its default.

export IPF_PROJECT_DIR="${IPF_PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
export IPF_METHSCAN_RUN_DIR="${IPF_METHSCAN_RUN_DIR:-${IPF_PROJECT_DIR}/Results/Methscan/CYL_ZCP_full_20260826_final}"
export IPF_INPUT_MANIFEST="${IPF_INPUT_MANIFEST:-${IPF_METHSCAN_RUN_DIR}/00_scanpy_selected/input_manifest.tsv}"
export IPF_FILTERED_CELL_IDS="${IPF_FILTERED_CELL_IDS:-${IPF_METHSCAN_RUN_DIR}/03_filtered/column_header.txt}"
export IPF_COV_DIR="${IPF_COV_DIR:-${IPF_PROJECT_DIR}/Data/30wcov}"
# Coverage and the historical ALLC directory are fallback inputs only. The
# maintained route resolves original indexed ALLCs from the MethSCAn manifest.
export IPF_EXISTING_ALLC_DIR="${IPF_EXISTING_ALLC_DIR:-/home/lijia/jiangyuanpei/methscan/xunyin/IPF_tissue/allcools_5kbin/input_allc}"
export IPF_ANNOTATION="${IPF_ANNOTATION:-${IPF_PROJECT_DIR}/Results/Scanpy/E_CYL_ZCP_notebook/cell_id_cell_type.tsv}"
export IPF_CHROM_SIZES="${IPF_CHROM_SIZES:-${IPF_PROJECT_DIR}/Supplementary/hg38.canonical.chrom.sizes}"
export IPF_BLACKLIST="${IPF_BLACKLIST:-${IPF_PROJECT_DIR}/Supplementary/ENCFF356LFX_GRCh38_blacklist.bed.gz}"
export IPF_BLACKLIST_MD5="${IPF_BLACKLIST_MD5:-393688b4f06c9ce26165d47433dd8c37}"
export IPF_BLACKLIST_FRACTION="${IPF_BLACKLIST_FRACTION:-0.2}"

# Keep the prior no-blacklist results intact. This root is required because the
# blacklist changes retained 5-kb features and therefore the MethylVI input.
export IPF_MVI_ROOT="${IPF_MVI_ROOT:-${IPF_PROJECT_DIR}/Results/MethylVI_30wcov_allcools_blacklist_f0p2}"
export IPF_ALLCOOLS_ROOT="${IPF_ALLCOOLS_ROOT:-${IPF_MVI_ROOT}/allcools_5kb}"
export IPF_ALLC_DIR="${IPF_ALLC_DIR:-${IPF_ALLCOOLS_ROOT}/input_allc}"
export IPF_ALLC_TABLE="${IPF_ALLC_TABLE:-${IPF_ALLCOOLS_ROOT}/selected_cells.allc.tsv}"
export IPF_MCDS="${IPF_MCDS:-${IPF_ALLCOOLS_ROOT}/mcg_5kb.mcds}"
export IPF_ALLCOOLS_H5AD="${IPF_ALLCOOLS_H5AD:-${IPF_ALLCOOLS_ROOT}/mcg_5kb.clustered.h5ad}"
export IPF_MVI_INPUT="${IPF_MVI_INPUT:-${IPF_MVI_ROOT}/ipf_allcools_5kb_methylvi_input.h5mu}"
export IPF_MVI_RESULTS="${IPF_MVI_RESULTS:-${IPF_MVI_ROOT}/results}"

# ALLCools and MethylVI use separate environments to avoid Python conflicts.
export IPF_ALLCOOLS_ENV="${IPF_ALLCOOLS_ENV:-/home/lijia/jiangyuanpei/miniforge3/envs/allcools}"
export IPF_METHYLVI_ENV="${IPF_METHYLVI_ENV:-/home/lijia/luozhixiong/miniconda3/envs/methylvi}"

export IPF_CELLTYPE_KEY="${IPF_CELLTYPE_KEY:-manual_celltype}"
# CYL/ZCP is the project-defined integration batch.
export IPF_BATCH_KEY="${IPF_BATCH_KEY:-sample_id}"

export IPF_BIN_SIZE="${IPF_BIN_SIZE:-5000}"
export IPF_MC_CONTEXT="${IPF_MC_CONTEXT:-CGN}"
export IPF_HYPO_SCORE_CUTOFF="${IPF_HYPO_SCORE_CUTOFF:-0.9}"
export IPF_BINARIZE_CUTOFF="${IPF_BINARIZE_CUTOFF:-0.95}"
# The 30k branch is built first; the nested 10k input is derived from its
# deterministic prevalence ranking without rescanning ALLCs.
export IPF_TARGET_FEATURES="${IPF_TARGET_FEATURES:-30000}"
export IPF_FEATURE_TARGETS="${IPF_FEATURE_TARGETS:-10000 30000}"
export IPF_LSI_COMPONENTS="${IPF_LSI_COMPONENTS:-100}"
export IPF_LSI_P_CUTOFF="${IPF_LSI_P_CUTOFF:-0.1}"
export IPF_ALLCOOLS_NEIGHBORS="${IPF_ALLCOOLS_NEIGHBORS:-25}"
export IPF_ALLCOOLS_LEIDEN_RESOLUTION="${IPF_ALLCOOLS_LEIDEN_RESOLUTION:-1.0}"
export IPF_CONSENSUS_LEIDEN_REPEATS="${IPF_CONSENSUS_LEIDEN_REPEATS:-500}"
export IPF_CONSENSUS_LEIDEN_RESOLUTION="${IPF_CONSENSUS_LEIDEN_RESOLUTION:-0.5}"

# yuanpei processes every coverage file. Cells without a manual label remain in
# the unsupervised analysis and receive manual_celltype=Unknown.
export IPF_EXPECTED_CELLS="${IPF_EXPECTED_CELLS:-6264}"
export IPF_INCLUDE_UNANNOTATED="${IPF_INCLUDE_UNANNOTATED:-1}"
export IPF_MAX_CELLS="${IPF_MAX_CELLS:-0}"
export IPF_THREADS="${IPF_THREADS:-16}"
export IPF_SEED="${IPF_SEED:-0}"
export IPF_EPOCHS="${IPF_EPOCHS:-500}"
export IPF_BATCH_SIZE="${IPF_BATCH_SIZE:-32}"

# Supervised UMAP is calculated from the fixed MethylVI latent representation;
# Unknown manual labels remain unsupervised.
export IPF_SUPERVISED_TARGET_KEY="${IPF_SUPERVISED_TARGET_KEY:-${IPF_CELLTYPE_KEY}}"
export IPF_SUPERVISED_TARGET_WEIGHTS="${IPF_SUPERVISED_TARGET_WEIGHTS:-0.2 0.5 0.7 0.9}"
export IPF_SUPERVISED_NEIGHBORS="${IPF_SUPERVISED_NEIGHBORS:-15}"
export IPF_SUPERVISED_MIN_DIST="${IPF_SUPERVISED_MIN_DIST:-0.5}"
