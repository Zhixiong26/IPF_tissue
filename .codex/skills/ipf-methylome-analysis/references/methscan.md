# MethSCAn workflow

Read this reference when implementing, checking, submitting, or interpreting the project MethSCAn stage.

## Starting contract

- Current formal source is the project-local `Data/ALLCools/` staging root, containing original CYL/ZCP `allcools.tar.gz` archives. Do not submit a run until both complete archives are present, readable, and their members/indices are validated on the intended compute node.
- The historic `/home/lijia/jiangyuanpei/methscan/xunyin/IPF_tissue/allcools_5kbin/input_allc` source is a 6,554-cell CGN-only `30wcov`-derived input. It is historical validation material only, not the current formal source and not evidence that the archive/RNA-gated workflow has passed.
- Use MethSCAn 1.1.0 at `/home/lijia/jiangyuanpei/miniforge3/envs/MethSCAn/bin/methscan`; do not modify that shared environment.
- Call `methscan prepare --input-format allc` directly. Canonical input links must yield `sample_barcode` as the MethSCAn cell name.
- Canonical RNA annotations come from `Results/Scanpy/E_CYL_ZCP_notebook/cell_id_cell_type.tsv`. Before any MethSCAn command, retain only ALLCs with an exact `cell_id` match whose RNA `cell_type` is neither empty nor literal `NA`. Write excluded cells and reasons to an audit table. This is the only cell-type selection: MethSCAn filter subsequently removes cells only by methylation QC, and VMR-Scanpy inherits labels from the entry manifest without post-filter cell-type rematching.

## Entry points

Use `Scripts/Methscan/01_prepare_allc_inputs.py` for discovery, archive extraction, validation, and candidate manifests; `01_select_scanpy_cells.py` for pre-prepare Scanpy whitelist selection; `02_prepare_methscan.py` for prepare; `run_methscan.sbatch` for the full chain; `03_vmr_scanpy.py` for VMR embedding/clustering; and `04_summarize_run.py` for final product validation.

The full chain is ALLC intake → Scanpy non-empty/non-`NA` whitelist → prepare → filter → smooth → scan → matrix → Scanpy PCA/neighbours/UMAP/Leiden → optional TSS profile → run summary. The default TSS input is the sort-only copy `Supplementary/human_hg38_TSS.methscan.bed` of the user-provided hg38 TSS BED; original coordinates and strand remain unchanged.

The project MethSCAn cell filter is `min-sites=300000`, `min-meth=50`, and `max-meth=100`. The site threshold is the technical covered-CpG eligibility rule. MethSCAn takes percentages and keeps equality at the minimum, so `min-meth=50` implements overall mCG ≥0.50; the 100% maximum is nonrestrictive. Do not add mCH/mCCC or a composite QC flag unless the user explicitly expands the contract. Keep code, README, Report, and the reference example synchronized.

## Execution decisions

- For the current contract, omit the external QC table and apply only the maintained MethSCAn site/overall-mCG thresholds. Do not imply that this also evaluates mapping, mCH, or mCCC.
- The former standalone per-cell QC workflow was removed and is not a default dependency. Design a new external QC-table gate only if the user explicitly expands the contract. Current CG-only ALLCs cannot supply mCH/mCCC.
- After archive staging and compute-node preflight pass, run a balanced small-cell smoke test before the full dataset. The same `IPF_METHSCAN_MAX_CELLS` applies to intake and Scanpy selection; record excluded smoke-subset cells.
- Require `stage_status.tsv`, per-stage portable resource JSON records, cell/VMR missingness tables, and strict equality/monotonicity checks across manifest, prepare, filter, matrix, and Scanpy cell counts before accepting `run_summary.json` as complete.
- Every run uses a new result root. Do not overwrite or delete partial MethSCAn data automatically; inspect failures and choose a new root unless a step is proven resumable.
- Verify the ALLC source from the intended compute node before submission. A visible login-node path is not sufficient.
- A submitted job is not a completed stage. Confirm `run_summary.json`, matrix products, H5AD, tables, figures, and logs before updating the stage to complete.
