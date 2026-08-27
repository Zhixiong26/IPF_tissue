# MethSCAn workflow

Read this reference when implementing, checking, submitting, or interpreting
the project MethSCAn stage.

## Current contract and entry points

- Source ALLCs are already extracted below `Data/ALLCools/`. Treat that tree
  as read-only: do not extract archives, copy inputs, or modify source ALLCs.
- The RNA gate is
  `Results/Scanpy/E_CYL_ZCP_notebook/cell_id_cell_type.tsv`. Retain an ALLC
  only when its canonical `<sample_id>_<17bp_barcode>` ID exactly matches an
  RNA `cell_id` whose `cell_type` is neither empty nor literal `NA`.
- The recommended production entry point is
  `Scripts/Methscan/submit_methscan_pipeline.sh`; run it once with one new
  directory below `Results/Methscan/`. It submits dependent Slurm jobs for
  selection/conversion, prepare, filter, smooth, three parallel threshold
  branches, and summary. The older `run_methscan.sbatch` remains a monolithic
  compatibility entry point.
- The runner calls, in order: `01_select_scanpy_cells.py`,
  `02_convert_allc_to_cov.py`, `03_prepare_methscan.py`, MethSCAn `filter` /
  `smooth` / `scan` / `matrix`, `04_vmr_scanpy.py`, and
  `05_summarize_run.py`. `06_run_with_resources.py` records elapsed time,
  return code, and resource usage for every stage.

The current implementation deliberately converts selected indexed ALLCs to
CpG-only Bismark `.cov.gz` files in the run's `01_cov/` directory, then calls
`methscan prepare --input-format bismark`. Do not replace this with a native
ALLC prepare invocation without an explicitly validated workflow change.

## Fixed defaults

`00_methscan_config.sh` is the source of truth for defaults. The current ones
are `min-sites=300000`, `min-meth=50`, and `max-meth=100`; this implements the
technical covered-CpG gate and overall mCG >= 0.50. It does not add an mCH,
mCCC, mapping-QC, or composite-QC gate.

After filtering, smooth with bandwidth 1,000 bp. Scan each variance threshold
`0.01`, `0.02`, and `0.05` with bandwidth 2,000 bp, step size 100 bp, and at
least six cells. For each branch, derive the residual matrix and run VMR
Scanpy with a 5% minimum region-cell fraction, at least 100 covered regions
per cell, 30 PCs, 15 neighbours, Leiden resolution 0.8, and seed 20260825.

The batch runner uses the verified MethSCAn and Scanpy interpreters, defaults
to 32 CPUs/256 GB on `fat`, and refuses submission unless the new result root
is below `Results/Methscan/`, inputs are available, and at least 500 GB is
free there. Scheduler stdout/stderr belong in `Scripts/Methscan/logs/`.

The split submission allocates 55 CPUs/250 GB on `cu03` to selection and
ALLC-to-COV, 4 CPUs/16 GB separately to prepare, filter, and smooth, and 18
CPUs/80 GB to each threshold branch. The three QC tools are serial; separating
their jobs avoids reserving the large conversion allocation for them.

## Run and completion rules

Use a fresh result root; never overwrite or resume a partial run by default.
Inspect `stage_status.tsv` and the stage resource JSON records after any
failure or cancellation.

Before accepting a run as complete, require `run_summary.json` with
`status=complete`, all three `VMRs.bed` files, all matrix products, each
branch's Scanpy H5AD/figures/tables, and successful ID checks. The summary
enforces that selected cells contain no empty/`NA` RNA labels; prepare IDs
match the selected manifest; filtered IDs are a subset of prepared IDs; matrix
columns equal filtered IDs; and each embedding is a subset of filtered IDs.
A submitted job or a populated partial directory is not completion evidence.

Use a bounded, balanced smoke run after changing this workflow, with a
separate temporary output root. Preserve source ALLCs and failed production
outputs for diagnosis.
