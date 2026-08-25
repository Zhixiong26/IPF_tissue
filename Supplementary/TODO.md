# Project TODO / 项目待办

## Pending / 待完成

- [ ] **恢复 full-context 逐细胞甲基化数据。** 当前 ALLC 来自六列 Bismark CpG coverage，只含 `CGN`，因此只能计算 mCG 和 covered CpG sites。需要从保留 CG、CH 和精确三核苷酸 context 的 methylation calls 重新生成逐细胞 ALLC，或提供等价的逐细胞 mCH/mCCC 统计。完成标准：`04_ALLC_methylation_QC.py` 对真实数据产生非缺失的 `mCH`、`mCCC`，并通过 `mc <= cov` 和 context 检查。

  **Recover full-context per-cell methylation data.** Current ALLCs originate from six-column Bismark CpG coverage and contain only `CGN`, so they support only mCG and covered CpG sites. Rebuild per-cell ALLCs from methylation calls retaining CG, CH, and exact trinucleotide contexts, or provide equivalent per-cell mCH/mCCC metrics. Completion criterion: `04_ALLC_methylation_QC.py` produces non-missing `mCH` and `mCCC` on real data and passes `mc <= cov` and context checks.
