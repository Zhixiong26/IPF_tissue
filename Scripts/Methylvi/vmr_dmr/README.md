# VMR + all unique pooled-DMR MethylVI

This workflow adds pooled CYL+ZCP cell-type DMR information to each existing
MethSCAn-VMR MethylVI input and retrains six models without modifying the
completed `MethylVI_clean6264_10k30k` results.

## Fixed analysis definition

- Source: all 91 pooled cell-type comparisons under
  `Methscan/CYL_ZCP_full_20260826_final/pooled/07_methdiff`.
- Retain every row with raw `p < 0.01` and `abs(meth_A - meth_B) >= 0.25`.
- Define the hypo cell type from MethSCAn's `low_group_label`.
- Remove exact duplicates within each hypo cell type and merge all overlapping
  retained DMR intervals. There is no Top-N/Top200 truncation.
- Apply the same hg38 blacklist rule used by the VMR workflow (remove source
  intervals with at least 20% blacklist overlap).
- Preserve each existing Top10k/Top30k VMR input. Append only merged DMRs that
  have zero genomic overlap with that model's selected VMRs, preventing the
  same CpG from being represented twice.
- Train `{0.01,0.02,0.05} x {Top10k,Top30k}` using integer `mc/cov`, the same
  6,264 cells, `sample_id` batch key, seed 0, and early stopping.

The completed preparation currently contains 1,104,139 exact-unique
hypo-DMRs and 321,910 merged non-overlapping DMR intervals. The merged regions
cover 967,355,000 bp and the maximum interval length is 58 kb.

## Submission

```bash
bash /home/lijia/luozhixiong/IPF_tissue/Scripts/Methylvi/vmr_dmr/submit_vmr_dmr_pipeline.sh
```

The submitter reuses `shared/dmr/prepare.COMPLETE`, builds the pooled-DMR count
matrix once, derives six model-specific non-overlap unions, trains the six
models, generates the standard cell-type/sample/condition/Leiden UMAPs, and
validates all outputs before writing `workflow.COMPLETE`.

All jobs are pinned to `cu03`. Count construction uses 40 CPU/48G; input joins
use 4 CPU/40G; training uses 32 CPU/48G. The six join/train pairs are serialized
to keep only one memory-intensive task active on the node at a time.

## Current production run

The DAG was submitted on 2026-09-17. DMR preparation is reused; count job
`309625` is running on cu03. Serialized join/train jobs are `309626–309637`,
and final validator `309638` waits for all six models. The workflow remains
`in progress` until `run_summary.json` and `workflow.COMPLETE` are created.

Results are written only below:

```text
Results/MethylVI_vmr_plus_all_unique_pooled_DMR/
```
