# MethSCAn workflow

Read this reference when implementing, checking, submitting, or interpreting the project MethSCAn stage.

## Starting contract

- Current verified source: `/home/lijia/jiangyuanpei/methscan/xunyin/IPF_tissue/allcools_5kbin/input_allc`. Generation logs show that these seven-column CGN-only ALLCs were converted from `30wcov`; preserve that provenance and do not call them original full-context ALLCs.
- This coverage-derived source is accepted for the current CpG-only VMR workflow because the requested methylation-level gate is overall mCG only. It is not suitable for mCH, mCCC, or strand-aware analyses. Original full-context ALLCs remain archived under `/mnt/data04/jiangyuanpei/xunyin_260727/data/allcools/<sample>/allcools.tar.gz` but are not required for the current contract.
- Current scope has 6,554 indexed ALLCs: CYL 3,165 and ZCP 3,389. Observed context is CGN only.
- Use MethSCAn 1.1.0 at `/home/lijia/jiangyuanpei/miniforge3/envs/MethSCAn/bin/methscan`; do not modify that shared environment.
- Call `methscan prepare --input-format allc` directly. Canonical input links must yield `sample_barcode` as the MethSCAn cell name.
- Canonical RNA annotations come from `Results/Scanpy/E_CYL_ZCP_notebook/cell_id_cell_type.tsv`. Before any MethSCAn command, retain only ALLCs with an exact `cell_id` match whose RNA `cell_type` is not literal `NA`. Write excluded cells and reasons to an audit table. After MethSCAn filter, require an exact subset match to that selected manifest; every retained cell must retain non-`NA` `rna_cell_type` fields.

## Entry points

Use `Scripts/Methscan/01_prepare_allc_inputs.py` for discovery, archive extraction, validation, and candidate manifests; `01_select_scanpy_cells.py` for pre-prepare Scanpy whitelist selection; `02_prepare_methscan.py` for prepare; `02_match_filtered_scanpy.py` for post-filter matching; `run_methscan.sbatch` for the full chain; `03_vmr_scanpy.py` for VMR embedding/clustering; and `04_summarize_run.py` for final product validation.

The full chain is ALLC intake → Scanpy non-`NA` whitelist → prepare → filter → post-filter Scanpy audit → smooth → scan → matrix → Scanpy PCA/neighbours/UMAP/Leiden → run summary. TSS profile is optional and stays disabled until the reference build and alphabetically sorted strand-aware BED are confirmed.

The project MethSCAn cell filter is `min-sites=300000`, `min-meth=50`, and `max-meth=100`. The site threshold is the technical covered-CpG eligibility rule. MethSCAn takes percentages and keeps equality at the minimum, so `min-meth=50` implements overall mCG ≥0.50; the 100% maximum is nonrestrictive. Do not add mCH/mCCC or a composite QC flag unless the user explicitly expands the contract. Keep code, README, Report, and the reference example synchronized.

## Execution decisions

- For the current contract, omit the external QC table and apply only the maintained MethSCAn site/overall-mCG thresholds. Do not imply that this also evaluates mapping, mCH, or mCCC.
- The former standalone per-cell QC workflow was removed and is not a default dependency. Design a new external QC-table gate only if the user explicitly expands the contract. Current CG-only ALLCs cannot supply mCH/mCCC.
- Always run a balanced small-cell smoke test before the full dataset. Measure peak memory, output size, VMR count, missingness, and retained cells before changing resource requests or thresholds.
- Require `stage_status.tsv`, per-stage portable resource JSON records, cell/VMR missingness tables, and strict equality/monotonicity checks across manifest, prepare, filter, matrix, and Scanpy cell counts before accepting `run_summary.json` as complete.
- Every run uses a new result root. Do not overwrite or delete partial MethSCAn data automatically; inspect failures and choose a new root unless a step is proven resumable.
- Verify the ALLC source from the intended compute node before submission. A visible login-node path is not sufficient.
- A submitted job is not a completed stage. Confirm `run_summary.json`, matrix products, H5AD, tables, figures, and logs before updating the stage to complete.
