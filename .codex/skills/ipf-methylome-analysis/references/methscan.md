# MethSCAn workflow

Read this reference when implementing, checking, submitting, or interpreting the project MethSCAn stage.

## Starting contract

- Start from per-cell ALLCools ALLC files or an archive containing them, not from the derived `Data/30wcov` coverage tree.
- Current verified source: `/home/lijia/jiangyuanpei/methscan/xunyin/IPF_tissue/allcools_5kbin/input_allc`.
- Current scope has 6,554 indexed ALLCs: CYL 3,165 and ZCP 3,389. Observed context is CGN only.
- Use MethSCAn 1.1.0 at `/home/lijia/jiangyuanpei/miniforge3/envs/MethSCAn/bin/methscan`; do not modify that shared environment.
- Call `methscan prepare --input-format allc` directly. Canonical input links must yield `sample_barcode` as the MethSCAn cell name.

## Entry points

Use `Scripts/Methscan/01_prepare_allc_inputs.py` for discovery, archive extraction, validation, optional QC selection, balanced smoke selection, and manifests. Use `02_prepare_methscan.py` for prepare, `run_methscan.sbatch` for the full chain, `03_vmr_scanpy.py` for VMR embedding/clustering, and `04_summarize_run.py` for final product validation.

The full chain is ALLC intake → prepare → filter → smooth → scan → matrix → Scanpy PCA/neighbours/UMAP/Leiden → run summary. TSS profile is optional and stays disabled until the reference build and alphabetically sorted strand-aware BED are confirmed.

## Execution decisions

- If per-cell QC is unfinished, omit the QC table and label the run accordingly; do not imply that all cells passed QC.
- When QC becomes available, review the selected Boolean column. Current CG-only ALLCs cannot supply mCH/mCCC, so a final flag depending on those metrics may be `NA` and must not be treated as pass.
- Always run a balanced small-cell smoke test before the full dataset. Measure peak memory, output size, VMR count, missingness, and retained cells before changing resource requests or thresholds.
- Every run uses a new result root. Do not overwrite or delete partial MethSCAn data automatically; inspect failures and choose a new root unless a step is proven resumable.
- Verify the ALLC source from the intended compute node before submission. A visible login-node path is not sufficient.
- A submitted job is not a completed stage. Confirm `run_summary.json`, matrix products, H5AD, tables, figures, and logs before updating the stage to complete.
