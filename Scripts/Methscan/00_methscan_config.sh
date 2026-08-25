#!/usr/bin/env bash
# Shared defaults for the preliminary CYL/ZCP MethSCAn workflow.

project_dir=/home/lijia/luozhixiong/IPF_tissue
methscan_exe=/home/lijia/jiangyuanpei/miniforge3/envs/MethSCAn/bin/methscan
methscan_python=/home/lijia/jiangyuanpei/miniforge3/envs/MethSCAn/bin/python
scanpy_python=/home/lijia/jiangyuanpei/miniforge3/envs/allcools/bin/python
rna_annotation_table=${IPF_METHSCAN_RNA_ANNOTATION_TABLE:-${project_dir}/Results/Scanpy/E_CYL_ZCP_notebook/cell_id_cell_type.tsv}
rna_exclude_cell_type=${IPF_METHSCAN_RNA_EXCLUDE_CELL_TYPE:-NA}

# Project-local ALLCools staging root. It contains one `allcools.tar.gz` archive
# per sample; intake safely expands these archives inside each new run directory.
allc_source=${IPF_METHSCAN_ALLC_SOURCE:-${project_dir}/Data/ALLCools}
samples=(CYL ZCP)

# Technical covered-CpG eligibility plus the requested overall mCG threshold.
# MethSCAn expresses methylation thresholds as percentages and treats the
# minimum as inclusive. max_meth=100 is only the valid-domain ceiling and does
# not impose an effective project upper filter.
min_sites=300000
min_meth=50
max_meth=100

smooth_bandwidth=1000
scan_bandwidth=2000
scan_stepsize=100
scan_var_threshold=0.02
scan_min_cells=6

threads=32
prepare_chunksize=10000000
min_free_gb=500
allc_validation_records=10000
allc_intake_workers=16
max_cells=${IPF_METHSCAN_MAX_CELLS:-0}

# Optional, only after the genome build and sorted BED have been confirmed.
tss_bed=${IPF_METHSCAN_TSS_BED:-}
tss_strand_column=6

# Downstream VMR representation defaults.
vmr_min_cell_fraction=0.05
cell_min_regions=100
pca_components=30
neighbors=15
leiden_resolution=0.8
random_seed=20260825
