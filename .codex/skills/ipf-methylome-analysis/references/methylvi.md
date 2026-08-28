# MethylVI routes

Read this reference before changing, running, or interpreting project MethylVI.

## Routes and entry points

| Route | Entry point | Upstream contract |
|---|---|---|
| ALLCools 5-kb | `Scripts/Methylvi/allcools/run.sh` | MethSCAn post-filter cell IDs mapped to original indexed ALLCs. |
| MethSCAn VMR | `Scripts/Methylvi/vmr/run.sh` | The same post-filter cells and one VMR branch from that run. |

The MethylVI input always has integer `mCG.layers['mc']` and `mCG.layers['cov']`.
Require `mc <= cov`, unique cell/feature IDs, and source/config manifests before
training. ALLCools hypo-scores are for region selection and clustering only.

Production uses `Scripts/Methylvi/submit_methylvi_8models.sh` once. It builds
ALLCools 5-kb × {10k,30k} and VMR {0.01,0.02,0.05} × {10k,30k}. All models use
the 6,264 IDs in `03_filtered/column_header.txt`; the selected manifest only
resolves original ALLC paths. Top 10k must be nested within top 30k. Train with
`sample_id` as batch key and require ordinary UMAP/Leiden, supervised UMAP,
sequencing-depth, overall-mCG, and mean-mCG diagnostics before completion.

## Current upstream location

The current formal MethSCAn run root is
`Results/Methscan/CYL_ZCP_full_20260826_final`. Its selected-ALLC manifest is
`00_scanpy_selected/input_manifest.tsv`; it records original indexed ALLCs
under `Data/ALLCools/`. The run is complete (`run_summary.json` has
`status=complete`) and provides VMR branches for `0.01`, `0.02`, and `0.05`.
The maintained eight-model experiment uses all three branches; verify every
`VMRs.bed` against this same run manifest before submission.

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

For future Slurm submissions, use approximately 2 GB RAM per requested CPU
unless measured peak RSS demonstrates a stage-specific exception. The current
maintained DAG defaults are builders 45 CPU/90G, trainers 45 CPU/90G, and final
summary 2 CPU/4G. Never assume edits to an sbatch file alter jobs that are
already queued or running.

The 2026-08-27 balanced 40-cell smoke passed 8/8 reduced models at nested
100/300 features and two epochs. AnnData `uns` mappings must not use arbitrary
cell-type labels as dictionary keys because labels can contain `/`, which is an
invalid HDF5 group key. Store such mappings as aligned label/code arrays; JSON
sidecars may retain the human-readable dictionary.

## Current formal run and plotting guidance

The formal output root is `Results/MethylVI_clean6264_10k30k`. The three VMR
builders completed and all six VMR models have `model.COMPLETE` markers. The
old trainer wrappers for jobs 307579/307580/307582/307583/307585/307586 exited
after training during plotting because the submitted snapshot lacked
`VMR_MVI_RESULTS`; recovery job 307591 successfully generated the ordinary
UMAP, supervised-UMAP, and methylation-QC outputs without retraining. Do not
interpret those historical exit codes as failed model fits; verify the model
marker, embedding, and plot summary instead.

ALLCools builder 307577 is the remaining bottleneck. It writes temporary Zarr
chunks before the final `mcg_5kb.mcds` assembly; the final file and completion
marker are the authoritative completion evidence. Jobs 307587/307588 (the
ALLCools 10k/30k models) and summary job 307592 are dependency-gated and must
not be marked complete while 307577 is running.

The standard supervised-UMAP weights are 0.2, 0.5, 0.7, and 0.9 for
sensitivity output. Because high target weights can dominate the latent-space
graph and create compressed/fragmented layouts, use 0.2 as the primary view;
0.3 or 0.4 may be generated from an existing embedding as auxiliary outputs,
without retraining. The six VMR model routes each contain 19 PNGs (ordinary
UMAP, QC, and four-weight supervised UMAP) after recovery.
