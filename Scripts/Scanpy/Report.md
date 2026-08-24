# Scanpy 分析报告 / Analysis report

状态 / Status：CYL/ZCP notebook 已准备，分析参数待交互确认 / CYL/ZCP notebook prepared; analysis parameters await interactive confirmation

CYL 和 ZCP 的本地 10x filtered matrix ZIP 已放入 `Data/Matrix/`：CYL 有 4,264 cells，ZCP 有 4,685 cells，均有 38,606 features。`E_CYL_ZCP_scanpy.ipynb` 以现有 `S12-2N.ipynb` 为模板，但不修改原 notebook，也不沿用 S12 的输入路径、cluster 编号或输出。

notebook 保留 `layers['counts']` 中的原始整数 counts，并在 `obs` 中保留 cohort 与 pre-filter QC。候选 QC 阈值、是否回归 covariates、PCA/邻接参数、Leiden resolution 和 cluster-to-cell-type mapping 均须在 notebook 中审核。保存 cell types/H5AD 的 cell 默认被 `ANALYSIS_CONFIRMED = False` 阻止。

现有 Python/Slurm 入口属于 notebook 验证前草稿，不作为正式分析提交。notebook 确认后再记录实际过滤数量、归一化/表示、邻接图、聚类、嵌入、QC 图、命令、环境及解释限制，并据此重写脚本。

No new run has been executed under `Results/Scanpy/`. When work begins, record input dimensions and metadata, filtering, normalization/representation, neighbours, clustering, embeddings, QC figures, commands, environments, and interpretation limits.
