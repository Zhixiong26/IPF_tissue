# ALLCools 5-kb → MethylVI report / 报告

## Maintained contract / 当前约定

The maintained cell list is the 6,264 IDs in the completed MethSCAn
`03_filtered/column_header.txt` (CYL 2,919; ZCP 3,345). The upstream selected
manifest is used only to resolve each ID to its original indexed ALLC. The
historical 6,554-cell `Data/30wcov` route is fallback code, not the production
cell-selection rule.

当前唯一细胞名单是 MethSCAn `03_filtered/column_header.txt` 中的 6,264 个
细胞（CYL 2,919；ZCP 3,345）。上游 manifest 只负责将 cell ID 映射到原始
indexed ALLC；历史 6,554-cell `Data/30wcov` 路线仅为兼容 fallback。

After 5-kb CGN MCDS generation and ENCODE GRCh38 blacklist filtering, bins are
ranked by the number of current cells passing the binarized hypo-score cutoff.
The production branch selects exactly top 30,000 bins, records the effective
`MVI_HYPO_PERCENT` and boundary ties in `feature_filter_summary.json`, and
derives a strictly nested top-10,000 H5MU using `selection_rank`.

5-kb CGN MCDS 与 blacklist 过滤后，按当前细胞集中 binarized hypo-score
阳性细胞数排序。正式流程精确选择 top 30,000，并记录动态
`MVI_HYPO_PERCENT`；top 10,000 从同一排序派生，不重复扫描 ALLC。

The production `Results/MethylVI_clean6264_10k30k` experiment completed both
ALLCools models (10k and 30k). Both have validated integer-count inputs,
`model.COMPLETE`, 20-dimensional embeddings, ordinary cell-type/sample/
condition/Leiden UMAPs, supervised UMAPs, and methylation QC outputs. They are
listed as `complete` in the experiment-level `run_summary.tsv`.

正式 `Results/MethylVI_clean6264_10k30k` 实验中的 ALLCools Top10k 和
Top30k 两个模型均已完成，并通过输入、模型标记、embedding、UMAP 与甲基化
QC 校验；实验级 `run_summary.tsv` 将两者均记录为 `complete`。
