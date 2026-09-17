# MethylVI report / MethylVI 报告

## Maintained workflow contract / 当前流程约定

All production routes use the same 6,264 MethSCAn-filtered cells (CYL 2,919;
ZCP 3,345), original indexed ALLCs, integer CGN methylated/coverage (`mc/cov`)
counts, and `sample_id` as the MethylVI batch key. Transformed methylation
scores and heatmap ratios are never used as MethylVI count matrices.

所有正式路线统一使用 MethSCAn filter 后的 6,264 个细胞（CYL 2,919；
ZCP 3,345）、原始 indexed ALLC、整数 CGN `mc/cov` 计数，并以
`sample_id` 作为 batch key。任何变换后的 score 或热图 ratio 都不会作为
MethylVI 训练计数。

The repository maintains three routes:

1. `allcools/`: ALLCools-selected 5-kb bins, Top10k/Top30k.
2. `vmr/`: MethSCAn VMR thresholds 0.01/0.02/0.05, each Top10k/Top30k.
3. `vmr_dmr/`: each selected VMR input plus all qualifying, merged,
   non-overlapping pooled cell-type DMRs.

## Completed 10k/30k baseline / 已完成基线

`Results/MethylVI_clean6264_10k30k` is complete (`workflow.COMPLETE`; final
summary job `307592`). Its experiment-level `run_summary.tsv` records 8/8
models as complete:

- ALLCools 5-kb: Top10k and Top30k;
- MethSCAn VMR: `{0.01,0.02,0.05} x {Top10k,Top30k}`.

Every model contains a saved MethylVI model, 20-dimensional latent embedding,
ordinary UMAPs colored by cell type/sample/condition/Leiden cluster,
supervised UMAPs, and methylation QC. The VMR source counts were 39,553,
80,818, and 166,618 before the MethylVI route's canonical/blacklist/coverage
selection.

## Train/validation split / 训练验证划分

The shared trainer does not define an independent test set. With scvi-tools
1.5.0.post1 defaults and seed 0, all 6,264 cells are globally shuffled (not
stratified by sample or cell type):

- training: 5,638 cells (CYL 2,638; ZCP 3,000);
- validation: 626 cells (CYL 281; ZCP 345);
- test: 0 cells.

Validation ELBO controls early stopping. Final latent embeddings and UMAPs are
computed for all 6,264 cells, so the UMAP is not an independent test-set
performance estimate.

## Pooled DMR upstream / pooled DMR 上游

MethSCAn pooled job `309188` completed on 2026-09-16. It combined CYL/ZCP cells
with the same RNA cell type and completed all 91 pairwise comparisons among 14
cell types using 6,264 cells. There were 0 failures and 0 fallback comparisons,
with 6,437,767 output DMR rows.

For the VMR+DMR route, every pooled DMR row satisfying raw `p < 0.01` and
`abs(methdiff) >= 0.25` is retained; there is no Top200 or other Top-N
truncation. After the shared 20% blacklist rule:

- 1,642,437 rows passed thresholds before exact deduplication;
- 1,104,139 DMRs were unique by hypo cell type and exact coordinates;
- overlap merging produced 321,910 non-overlapping DMR intervals;
- merged span: 967,355,000 bp; maximum interval length: 58 kb.

For each model, the selected Top10k/Top30k VMR intervals remain unchanged and
only pooled DMR intervals with zero genomic overlap to those VMRs are appended,
avoiding duplicate representation of the same CpGs.

## VMR+DMR run status — 2026-09-17 / 当前运行状态

Results are isolated under
`Results/MethylVI_vmr_plus_all_unique_pooled_DMR`; the completed baseline is
not overwritten. The workflow is currently **in progress** on `cu03`:

- `309625`: build pooled-DMR integer counts, running with 40 CPU/48G;
- `309626/309627`: var0.01 Top10k join/train;
- `309628/309629`: var0.01 Top30k join/train;
- `309630/309631`: var0.02 Top10k join/train;
- `309632/309633`: var0.02 Top30k join/train;
- `309634/309635`: var0.05 Top10k join/train;
- `309636/309637`: var0.05 Top30k join/train;
- `309638`: final validation and summary.

The six join/train pairs are serialized. Join jobs request 4 CPU/40G; training
jobs request 32 CPU/48G; summary requests 2 CPU/8G. Each trainer automatically
exports standard UMAPs for `cell_type`, `sample_id`, `condition`, and
`methylVI_leiden`, followed by supervised UMAP and methylation QC. Only final
job `309638` may create `workflow.COMPLETE`; pending jobs must not be described
as complete before that marker and `run_summary.json` exist.
