# CYL/ZCP Scanpy Notebook 操作契约

本文件面向执行该分析的人或自动化 skill，描述长期稳定的目标、边界和验证规则。运行状态、失败证据和修正记录写入 `Report.md`，不要在两处重复维护易变化的参数。

## 适用范围

唯一分析入口：

- `Notebooks/E_CYL_ZCP_scanpy.ipynb`

该 Notebook 负责读取 CYL/ZCP 10x filtered matrices，完成 QC、Scrublet、归一化、HVG、PCA、Harmony、邻接图、UMAP、Leiden、marker、人工注释和结果导出。

以下内容不属于当前工作流：独立 Python 脚本、Slurm wrapper、历史脚本结果和其他项目的 cluster 编号映射。

## 不变量与权限边界

- 输入位于项目 `Data/Matrix/`，分析期间只读。
- 正式输出位于 `Results/Scanpy/E_CYL_ZCP_notebook/`。
- 原始整数表达应保留在 `layers['counts']`。
- 不把其他样本或历史运行的 cluster 编号直接映射到当前细胞。
- 用户已授权在本工作流的自迭代阶段修改 canonical Notebook 中与 PCA、Harmony、邻接图、UMAP 和 Leiden 有关的参数并重跑。该授权不包含静默改变输入、QC 标准、marker 生物学判定或覆盖正式数据。
- 未经人工复核，不把探索性运行标记为正式结果，也不覆盖已有 H5AD、逐细胞注释或参数文件。
- 不通过放宽验证、删除失败证据或追求历史数字来制造“通过”。

## 文件角色

- `README.md`：稳定操作契约，仅在工作流目标或约束改变时更新。
- `Report.md`：当前状态、基线、每轮验证和修正决策；每轮工作后更新。
- Notebook：分析实现和参数的唯一事实来源。
- `confirmed_parameters.json`：一次正式运行实际采用的机器可读参数。
- `figure_manifest.json`：一次运行实际导出的图片清单。

## Jupyter Notebook 具体分析流程

以下流程对应 `Notebooks/E_CYL_ZCP_scanpy.ipynb` 当前实现。表中的参数是当前值；实际执行时仍以 Notebook 单元和本轮生成的 `confirmed_parameters.json` 为准。

### 1. 初始化与输入读取

- 导入 Scanpy、AnnData、Scrublet、Harmony、Pandas、NumPy 和绘图库，固定随机种子 `RANDOM_SEED = 0`。
- 从 `Data/Matrix/` 分别读取 CYL、ZCP 的 10x filtered matrix ZIP。
- 在临时目录解压 ZIP，并使用 `sc.read_10x_mtx(..., var_names='gene_symbols', make_unique=True)` 读取标准 `filtered_feature_bc_matrix/`。
- 给每个 barcode 添加 `CYL_` 或 `ZCP_` 前缀，在 `obs['cohort']` 记录来源，并检查输入路径和 cell ID。
- 将未归一化整数矩阵复制到 `layers['counts']`；按基因名标记线粒体基因 `MT-` 和核糖体基因 `RPL/RPS`。

### 2. 过滤前 QC 探索

- 使用 `sc.pp.calculate_qc_metrics()` 计算 `total_counts`、`n_genes_by_counts`、`pct_counts_mt` 和核糖体比例等指标。
- 每个 cohort 独立绘制最高表达的 20 个基因、QC 小提琴图、`total_counts–pct_counts_mt` 和 `total_counts–n_genes_by_counts` 散点图。
- 根据两个 cohort 的实际分布审核阈值，不从其他 Notebook 自动继承阈值。

### 3. 基础 QC 与基因过滤

每个 cohort 独立应用当前候选阈值：

| 参数 | 当前值 | 判定 |
|---|---:|---|
| `MIN_GENES` | 200 | `n_genes_by_counts >= 200` |
| `MAX_GENES` | 6000 | `n_genes_by_counts < 6000` |
| `MIN_COUNTS` | 500 | `total_counts >= 500` |
| `MAX_MT_PERCENT` | 5.0 | `pct_counts_mt < 5.0` |
| `MIN_CELLS_PER_GENE` | 3 | cohort 内 `filter_genes(min_cells=3)` |

Notebook 把基础 QC 判定写入 `obs['pass_basic_qc']`，记录各 cohort 过滤前后细胞数。若任一 cohort 在基础 QC 后少于 100 个细胞，则停止，避免不稳定的 doublet 推断。

### 4. 分 cohort Scrublet

- Scrublet 只使用通过基础 QC 的原始 `layers['counts']`，在合并样本前分别运行。
- 预期 doublet 率按 `0.004 × 当前 cohort 细胞数 / 1000` 动态计算。
- 当前使用 `n_prin_comps=30`、`use_approx_neighbors=False` 和随机种子 `0`。
- 将连续分数和自动判定分别写入 `obs['doublet_score']`、`obs['predicted_doublet']`，保存分数直方图并人工检查阈值。
- 排除 `predicted_doublet=True` 的细胞，只将 singlets 送入下游。

### 5. 合并、归一化与高变基因

- 使用 `anndata.concat(..., join='outer', merge='same', index_unique=None)` 合并各 cohort singlets，并再次确认 cell ID 唯一。
- 使用 `normalize_total(target_sum=10000)` 和 `log1p()` 归一化。
- 将全基因 log-normalized 表达保存到 `adata.raw`，供 marker 检验和绘图使用。
- 使用 `highly_variable_genes(n_top_genes=2000, batch_key='cohort', flavor='seurat')` 选择 HVG，并保存 HVG 诊断图。

### 6. Scale、PCA 与 Harmony

- 下游工作对象只保留 2,000 个 HVG。
- 当前 `REGRESS_COVARIATES=False`；只有明确审核后才回归 `total_counts` 和 `pct_counts_mt`。
- 使用 `scale(max_value=10)` 标准化，再以 ARPACK、50 个主成分和随机种子 `0` 计算 PCA。
- 保存 PCA 方差图和 Harmony 前按 cohort 着色的 PCA 图。
- 在 `X_pca` 上按 `cohort` 执行 Harmony，结果写入 `obsm['X_pca_harmony']`；当前最大迭代数为 20，初始 cluster 数按细胞数动态计算，随机种子为 `0`。

### 7. 邻接图、UMAP 与 Leiden

- 使用 `N_PCS=30`、`N_NEIGHBORS=15` 分别从原始 `X_pca` 和 Harmony 后 `X_pca_harmony` 构建邻接图。
- 分别保存 `X_umap_before_harmony` 和 `X_umap_after_harmony`，用相同参数比较 Harmony 前后 sample 分布。
- 在 Harmony UMAP 上检查 `total_counts`、`pct_counts_mt` 是否仍主导嵌入结构。
- 在 Harmony 后邻接图上使用 `LEIDEN_RESOLUTION=0.8`、随机种子 `0` 运行 Leiden，并检查每个 cluster 的 CYL/ZCP 构成。

### 8. Marker 检验与人工注释

- 使用 `rank_genes_groups(groupby='leiden', method='wilcoxon', use_raw=True)` 计算每个 cluster 的 marker，并绘制 Top 20 marker。
- 结合 ranked markers、经典肺细胞 marker、UMAP 位置、样本构成、QC 和稀有群证据审核注释。
- `cluster_to_cell_type` 仅适用于当前已经审核的 cluster 集合。实际 cluster 如有新增、缺失或重编号，Notebook 会停止，禁止静默复用旧映射。
- 将 Leiden cluster 映射到 `obs['cell_type']`，同一一级 cell type 可以合并多个 Leiden cluster，但原 cluster 身份仍保留。

### 9. 注释图与专项复核

- 从 `adata.raw` 绘制全局 cell-type marker dotplot；点大小表示表达比例，颜色表示平均表达。
- 单独绘制上皮相关 cluster 和稀有 cluster 的 marker dotplot，并汇总稀有群的细胞数、QC、doublet score 和 cohort 构成。
- 在 Harmony 后坐标分别生成 Leiden、cell type、sample 和 QC UMAP；只有存在经过验证的 `obs['group']` 时才生成 group UMAP。
- 所有 Notebook 图片以 300 dpi PNG 写入 `Results/Scanpy/E_CYL_ZCP_notebook/figures/`，同时登记到 manifest。

### 10. 确认、保护与导出

- 最终保存单元首先检查 `ANALYSIS_CONFIRMED`；该值只能在阈值、降维、聚类、marker 和注释均已审核后启用。
- 当前 `OVERWRITE_DATA_OUTPUTS=False`。正式输出已存在时，Notebook 会比较 cell ID 和逐细胞元数据；不一致则停止，不自动覆盖。
- 正式输出包括：

  - `rna_e_cyl_zcp_confirmed.h5ad`
  - `cell_id_cell_type.tsv`
  - `confirmed_parameters.json`
  - `figure_manifest.json`
  - `figures/`

- H5AD、逐细胞表、参数 JSON、图片清单和实际图片必须通过 README 后文规定的输出一致性检查。

## 执行闭环

### 1. 建立基线

执行前记录：

- Notebook SHA-256 和 Git 状态。
- 输入路径、存在性、大小及可获得的校验信息。
- 输出目录中已存在的正式产物。
- Notebook 当前的 `ANALYSIS_CONFIRMED` 与 `OVERWRITE_DATA_OUTPUTS` 值。
- Python kernel 和关键包版本。

若工作区已有修改，将其视为用户工作。不得重置、覆盖或混入本轮修正。

### 2. 干净重放

在兼容的 Jupyter kernel 中从第一单元顺序执行到最后一单元。诊断性重放优先保存为新的执行副本或新的结果位置，canonical Notebook 保持不变。

先验证上游分析，再决定是否允许最终保存。不要因为最终保存保护主动停止，就把已成功完成的上游分析误判为计算失败。

### 3. 分层验证

按以下顺序检查，上一层失败时不宣称下一层通过：

1. **执行完整性**：所有预期单元已运行；无未解释的 traceback；kernel 与依赖可识别。
2. **结构完整性**：cell ID 唯一；必要的 `obs`、`layers`、`obsm`、`uns` 字段存在；矩阵与元数据维度一致。
3. **分析合理性**：QC 流转可解释；Harmony 前后图、Leiden、marker、样本构成和最终注释 UMAP 经过复核。
4. **注释完整性**：每个实际 cluster 都有当前 marker 证据；新增、缺失或重编号 cluster 必须重新审核；最终注释 UMAP 与 marker、邻域和样本证据不存在明显冲突。
5. **输出一致性**：H5AD 细胞数与逐细胞表一致；参数 JSON 对应本轮设置；manifest 数量和文件实际存在性一致。

历史基线只能用于发现差异，不能替代本轮验证。

### 4. 最终注释 UMAP 与参数自迭代

完成初步注释后，必须检查最终 `cell_type` UMAP，而不是在代码无报错或标签表完整时结束。至少检查：

- 同一 cell type 是否出现无法由已知状态解释的严重碎裂。
- marker 明显不相容的 cell types 是否大面积混合，或存在异常桥接结构。
- 每个 cell type 内是否仍主要按 sample/cohort 分离，且这种分离是否有已知生物学依据。
- 稀有群是否稳定、是否被大群吞并，以及是否具有独立 marker 支持。
- UMAP 结构是否主要由 `total_counts`、`pct_counts_mt`、doublet score 或其他技术指标驱动。
- Leiden 边界、cell-type 合并结果、marker dotplot 与局部邻域是否相互一致。

UMAP 是二维可视化，不单独构成注释正确性的证据。不能为了获得更分离、更圆或更美观的图而牺牲 marker 特异性、合理的连续状态、稀有群或样本生物学差异。

若发现明确优化空间，可自主形成小规模、有界的候选参数，并按一次只改变一个参数家族的原则迭代：

| 参数家族 | 可调整内容 | 必须重跑的范围 |
|---|---|---|
| 特征与 PCA | HVG 数、是否回归协变量、PCA 成分数 | 从归一化/HVG 或 PCA 开始，重跑全部下游 |
| Harmony | batch key 已确认前提下的 `nclust`、`sigma`、迭代参数 | Harmony、邻接图、UMAP、Leiden、marker、注释 |
| 邻接图 | `N_PCS`、`N_NEIGHBORS` | 邻接图、UMAP、Leiden、marker、注释 |
| UMAP | `min_dist`、`spread` 和其他布局参数 | UMAP 与图形复核；不能据此声称 cluster 或注释已改善 |
| Leiden | `resolution` | Leiden、marker、cluster 映射和全部注释复核 |

每轮迭代遵循：

1. 保存当前基线的参数、图、cluster 数、marker 摘要、样本构成和输出路径。
2. 说明要改善的具体缺陷及候选参数的理由。
3. 将候选运行写入独立的迭代目录，不覆盖当前正式结果。
4. 从最早受影响步骤开始重跑；cluster 集合变化时废弃旧编号映射并重新计算 marker、重新注释。
5. 使用固定配色、类别顺序和绘图尺寸并排比较候选与基线。
6. 只有目标缺陷改善，且 marker、QC、样本结构、稀有群和稳定性没有明显退化时，才接受候选为新基线。
7. 更新 `Report.md`，继续检查新基线；仍有证据支持的优化空间时进入下一轮。

参数优化在以下情况视为收敛：最终注释 UMAP 没有可由参数修正的明显结构问题；marker 与邻域证据一致；连续两轮候选均未带来实质改善；或继续调整只改变视觉布局而不改善分析证据。

### 5. 诊断与窄修正

失败时先保留原始错误和最小复现证据，再归类：

- 环境/依赖：kernel、包版本、资源或文件格式问题。
- 输入：路径、内容、样本身份或矩阵结构变化。
- 执行状态：乱序执行、残留变量、部分输出或旧缓存。
- 分析逻辑：参数、数据变换、随机性或 API 行为变化。
- 生物学审核：cluster 结构、marker 或注释证据变化。

每次只修正有直接证据支持的最小范围，然后从受影响步骤之前重新执行，并复查全部下游不变量。PCA、Harmony、邻接图、UMAP 和 Leiden 参数可按上一节授权自主迭代；若修正需要改变其他 canonical Notebook 逻辑、覆盖正式数据或改变缺少证据支持的生物学判定，先停止并请求用户授权。

### 6. 记录与停止条件

每轮结束都在 `Report.md` 记录证据、结论和下一步。满足以下任一条件即停止自动迭代并交由用户决定：

- 需要修改已授权参数家族以外的 canonical Notebook 逻辑。
- 需要覆盖已确认的数据产品。
- 输入身份或预期分析目标不明确。
- cluster/marker 差异无法从表达证据得到可靠生物学判断。
- 连续两轮参数候选均未带来实质改善，或只改变视觉外观。
- 继续运行会显著增加资源消耗或产生不可逆影响。

## 正式完成标准

只有同时满足以下条件，才可在报告中写为“通过”：

- 执行、结构、分析、注释和输出五层验证均有证据，最终注释 UMAP 已完成合理性审查和必要的参数迭代。
- 所有实际 cluster 已基于本轮结果审核。
- 保存行为符合 `ANALYSIS_CONFIRMED` 和 `OVERWRITE_DATA_OUTPUTS` 的明确意图。
- 机器可读输出与 Notebook 本轮状态一致。
- `Report.md` 已记录本轮输入、代码身份、结果、差异和未决风险。
