# VMR + all unique pooled-DMR MethylVI report

## Input definition

- Cells: 6,264 MethSCAn-filtered cells (CYL 2,919; ZCP 3,345).
- Pooled DMR source: 91/91 completed comparisons from MethSCAn job `309188`.
- DMR filter: raw `p < 0.01` and `abs(methdiff) >= 0.25`.
- No Top200 or other Top-N truncation.
- Exact deduplication key: hypo cell type + chromosome + start + end.
- Blacklist: remove source intervals with at least 20% overlap.
- DMR merge: merge all overlapping retained intervals.
- VMR/DMR overlap: retain selected VMRs and append only zero-overlap DMRs.

Preparation completed with 1,642,437 qualifying rows, 1,104,139 exact-unique
hypo-DMRs, and 321,910 merged intervals spanning 967,355,000 bp. The longest
merged interval is 58 kb.

## Models and outputs

Six models are defined as MethSCAn VMR
`{0.01,0.02,0.05} x {Top10k,Top30k}`, with the model-specific non-overlapping
pooled DMR set appended. Each model uses integer CGN `mc/cov`, batch key
`sample_id`, seed 0, up to 500 epochs, and validation-ELBO early stopping.

Each completed route must contain:

- `input.h5mu` and `input.summary.json`;
- saved MethylVI model and `methylvi_embedding.h5ad`;
- cell-type, sample, condition, and Leiden UMAP PNGs;
- supervised UMAP and methylation QC outputs;
- `input.COMPLETE` and `model.COMPLETE`.

## Current status — 2026-09-17

The workflow is running serially on cu03:

| Stage | Job(s) | Resources | Status |
|---|---|---|---|
| Pooled-DMR counts | `309625` | 40 CPU / 48G | running |
| Six input joins | `309626, 309628, 309630, 309632, 309634, 309636` | 4 CPU / 40G each | dependency |
| Six MethylVI trains/plots | `309627, 309629, 309631, 309633, 309635, 309637` | 32 CPU / 48G each | dependency |
| Final validation | `309638` | 2 CPU / 8G | dependency |

Result root:
`Results/MethylVI_vmr_plus_all_unique_pooled_DMR`.

Status is **in progress**. Only job `309638` may write the final
`run_summary.json` and `workflow.COMPLETE` marker after all six models pass.
