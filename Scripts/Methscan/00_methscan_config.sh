#!/usr/bin/env bash
# Shared defaults for the preliminary CYL/ZCP MethSCAn workflow.

project_dir=/home/lijia/luozhixiong/IPF_tissue
methscan_exe=/home/lijia/jiangyuanpei/miniforge3/envs/MethSCAn/bin/methscan
methscan_python=/home/lijia/jiangyuanpei/miniforge3/envs/MethSCAn/bin/python
scanpy_python=/home/lijia/jiangyuanpei/miniforge3/envs/allcools/bin/python

# Starting point: a directory of ALLCools ALLC files, or an allcools tar archive.
allc_source=${IPF_METHSCAN_ALLC_SOURCE:-/home/lijia/jiangyuanpei/methscan/xunyin/IPF_tissue/allcools_5kbin/input_allc}
samples=(CYL ZCP)

# These MethSCAn-level filters are provisional defaults inherited from
# Scripts/Methscan/Methscan.md. Project QC selection can be supplied separately.
min_sites=300000
min_meth=20
max_meth=85

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
