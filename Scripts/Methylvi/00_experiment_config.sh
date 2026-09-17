#!/usr/bin/env bash
# Shared configuration for the clean-6264-cell, 10k/30k MethylVI experiment.
export MVI_PROJECT_DIR="${MVI_PROJECT_DIR:-/home/lijia/luozhixiong/IPF_tissue}"
export MVI_EXPERIMENT_ROOT="${MVI_EXPERIMENT_ROOT:-${MVI_PROJECT_DIR}/Results/MethylVI_clean6264_10k30k}"
export MVI_FEATURE_TARGETS="${MVI_FEATURE_TARGETS:-10000 30000}"
export MVI_VMR_THRESHOLDS="${MVI_VMR_THRESHOLDS:-0.01 0.02 0.05}"
export MVI_MAX_FEATURES="${MVI_MAX_FEATURES:-30000}"
export MVI_METHYLVI_ENV="${MVI_METHYLVI_ENV:-/home/lijia/luozhixiong/miniconda3/envs/methylvi}"
export MVI_ALLCOOLS_ENV="${MVI_ALLCOOLS_ENV:-/home/lijia/jiangyuanpei/miniforge3/envs/allcools}"
