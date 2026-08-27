# MethylVI routes

Read this reference before changing, running, or interpreting project MethylVI.

## Routes and entry points

| Route | Entry point | Upstream contract |
|---|---|---|
| ALLCools 5-kb | `Scripts/Methylvi/allcools/run.sh` | Historical six-column coverage and indexed CGN ALLCs. |
| MethSCAn VMR | `Scripts/Methylvi/vmr/run.sh` | One completed MethSCAn run's selected-ALLC manifest and one VMR branch from that same run. |

The MethylVI input always has integer `mCG.layers['mc']` and `mCG.layers['cov']`.
Require `mc <= cov`, unique cell/feature IDs, and source/config manifests before
training. ALLCools hypo-scores are for region selection and clustering only.

## Current upstream location

The current formal MethSCAn run root is
`Results/Methscan/CYL_ZCP_full_20260826_final`. Its selected-ALLC manifest is
`00_scanpy_selected/input_manifest.tsv`; it records original indexed ALLCs
under `Data/ALLCools/`. Its scan/VMR stages are not complete yet. Therefore do
not submit VMR-MethylVI until both `run_summary.json` says `status=complete`
and the chosen `04_scan/var_<threshold>/VMRs.bed` exists.

For another completed run, set `VMR_METHSCAN_RUN_DIR` and
`VMR_METHSCAN_VARIANCE`, or set `IPF_METHSCAN_VMR_SOURCE` and
`VMR_INPUT_MANIFEST` together. Never mix a VMR BED and selected-ALLC manifest
from different runs.

## Safe operation

Use `verify` before each stage and a bounded smoke test before a new full run.
Use a fresh result root when cell selection, regions, or parameters change.
Scheduler logs belong in `Scripts/Methylvi/logs/`. A submitted scheduler job is
not a completed analysis; verify the output manifest, embedding, and QC before
recording completion in the stage report.
