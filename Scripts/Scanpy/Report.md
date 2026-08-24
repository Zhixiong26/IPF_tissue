# Scanpy 分析报告 / Analysis report

状态 / Status：E 转录组 Scanpy 流程已准备，尚未正式运行 / transcriptome E Scanpy workflow prepared, not yet run

CYL 和 ZCP 的本地 10x filtered matrix ZIP 已放入 `Data/Matrix/`：CYL 有 4,264 cells，ZCP 有 4,685 cells，均有 38,606 features。`run_e_scanpy.sbatch` 在 `fat` partition 使用 16 CPUs/64 GiB，输出限于 `Data/Scanpy/`。

默认流程保留 `layers['counts']` 中的原始整数 counts，并在 `obs` 中保留 cohort 与 pre-filter QC。它不会把甲基化的 `manual_celltype` 自动套用到 RNA barcodes；RNA 注释将基于聚类、marker genes 和后续人工审核。

开始后在此记录实际过滤数量、归一化/表示、邻接图、聚类、嵌入、QC 图、命令、环境及解释限制。

No new run has been executed under `Results/Scanpy/`. When work begins, record input dimensions and metadata, filtering, normalization/representation, neighbours, clustering, embeddings, QC figures, commands, environments, and interpretation limits.
