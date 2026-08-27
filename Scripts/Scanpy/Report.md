# CYL/ZCP Scanpy Notebook 迭代报告

本文件是可持续更新的运行账本。最新状态放在最前；历史记录保留结论和证据，不复制 Notebook 中的完整实现。稳定规则见 `README.md`。

## 最新状态

- 状态：四项优先改进已实现；baseline 与隔离 candidate smoke test 均已无错误完成。
- Canonical Notebook：`Notebooks/E_CYL_ZCP_scanpy.ipynb`
- 本轮 Notebook 行为：已修改并完整执行；17 个代码单元顺序执行且无 error output。完整执行副本已保存到 Results，canonical Notebook 随后清为零输出交付状态。
- 正式输出目录：`Results/Scanpy/E_CYL_ZCP_notebook/`
- 已知保存开关：Notebook 当前文本中 `ANALYSIS_CONFIRMED = True`、`OVERWRITE_DATA_OUTPUTS = False`；每次运行前必须重新读取确认。
- 已授权范围：可修改并重跑 PCA、Harmony、邻接图、UMAP、Leiden 参数；候选运行不得覆盖正式结果。
- 当前未决事项：暂无实现阻塞；真实参数优化时应使用新的 `ITERATION_ID`，并比较候选 UMAP、marker 和审计表。

## 已知参考基线

以下是既有 Notebook 重放记录，用于差异检测，不是强制成功条件：

| 项目 | 已记录值 | 使用方式 |
|---|---:|---|
| CYL 输入 cells | 4,264 | 输入一致时比较 |
| ZCP 输入 cells | 4,685 | 输入一致时比较 |
| QC 与 doublet 后 cells | 8,696 | 差异出现时定位 QC/Scrublet 阶段 |
| Leiden clusters | 18（0–17） | 仅作结构变化警报；不得据此强套标签 |
| 随机种子 | 0 | 与实际参数 JSON 交叉检查 |

已记录的正式产物：

- `Results/Scanpy/E_CYL_ZCP_notebook/rna_e_cyl_zcp_confirmed.h5ad`
- `Results/Scanpy/E_CYL_ZCP_notebook/cell_id_cell_type.tsv`
- `Results/Scanpy/E_CYL_ZCP_notebook/confirmed_parameters.json`
- `Results/Scanpy/E_CYL_ZCP_notebook/figure_manifest.json`
- `Results/Scanpy/E_CYL_ZCP_notebook/figures/`

这些路径表示预期接口，不证明文件在下一次运行时仍存在或仍有效。

## 当前验证矩阵

| 层级 | 状态 | 当前证据 | 下一步 |
|---|---|---|---|
| 执行完整性 | 通过 | baseline 与 candidate 的 17 个代码单元均按 1–17 执行，无 error output | 真实候选继续同样检查 |
| 结构完整性 | 通过 | baseline 与 candidate 均为 8,696 cells、18 clusters；candidate H5AD/TSV 结构一致 | 参数改变后重新检查 |
| 分析合理性 | 通过 | baseline 重放生成 20 张图；核心细胞数和 cluster 结构保持一致 | 后续候选需单独比较 UMAP/marker |
| 注释完整性 | 通过 | baseline 为 `reviewed`；candidate 生成 18 行审计，17 个 proposed、1 个 ambiguous/Unassigned，状态为 `candidate_requires_review` | 候选不得直接升级为正式标签 |
| 输出一致性 | 通过 | baseline manifest 20 张图；candidate manifest 18 张图且零缺失；正式三个数据文件哈希在 smoke test 前后不变 | 每个真实候选使用新目录 |

状态只使用：`未验证`、`通过`、`失败`、`阻塞`。没有证据时不得写为“通过”。

## 迭代记录

### 2026-08-27：实现候选参数与注释自迭代后端

- 触发：终检发现错误输出污染、固定输出目录、UMAP 参数未集中和 cluster 变化后缺少候选注释流程。
- 修改：新增 baseline/candidate 运行模式、唯一迭代目录、集中式 PCA/Harmony/neighbors/UMAP/Leiden 参数、marker-score 候选注释与审计输出。
- Notebook 清洁：清除 16 个 `notebook controller is DISPOSED` 输出，并用 `ipf-allcools` kernel 完整重放。
- baseline 验证：17 个代码单元顺序完成、无错误；三个正式数据文件通过逐细胞保护而未覆盖；20 张图重新生成。
- 执行证据：完整副本保存为 `Results/Scanpy/E_CYL_ZCP_notebook/E_CYL_ZCP_scanpy_executed.ipynb`；canonical Notebook 验证后再次清除输出，避免提交内嵌图和 kernel 状态。
- 结果身份：manifest 写入 `run_kind=baseline`、`iteration_id=null`、`annotation_status=reviewed`。
- candidate smoke test：使用相同分析参数但强制进入 candidate 分支，独立写入 `/tmp/ipf_scanpy_candidate_smoke/iterations/same_params_smoke/`。
- candidate 结果：17 个 cluster 得到 marker-score proposal；cluster 13 因 Mast/T top1-top2 margin 仅 0.058，小于 0.20 阈值而保守标为 `Unassigned`。
- 隔离验证：candidate 使用 `candidate_cell_type`，没有正式 `cell_type` 列；H5AD/TSV/JSON 和 18 张图均使用候选文件名，正式数据哈希未变化。
- 结论：baseline 保护、候选隔离、候选审计与低置信拒绝机制均工作正常。

### 2026-08-27：扩展注释 UMAP 自迭代范围

- 触发：用户要求迭代不能只处理报错和重复确认注释，还要评价最终注释 UMAP，并在存在明显优化空间时自主调参、重跑和再次判断。
- 授权：允许修改 PCA、Harmony、邻接图、UMAP 和 Leiden 参数；不包含静默改变输入、QC、生物学结论或覆盖正式结果。
- 修改：README 新增 UMAP 合理性检查、可调参数家族、受影响重跑范围、候选比较、接受标准和收敛条件。
- 验证：本轮仅更新文档，没有执行 Notebook 或产生候选分析。
- 结论：后续 skill 应持续迭代到分析证据收敛，而不是以无报错或已有 cell-type 标签作为终点。

### 2026-08-27：文档结构修正

- 触发：需要让 README/Report 可被 skill 稳定使用，并支持自我迭代和纠错。
- 观察：旧文档能说明入口和输出，但没有失败分类、修正边界、停止条件或统一记录格式。
- 修改：README 改为稳定操作契约；Report 改为状态账本，增加参考基线、验证矩阵、修正记录和迭代模板。
- 验证：核对 Notebook 中的输出目录、随机种子、保存开关和正式输出文件名；未修改 Notebook。
- 结论：文档层面的缺口已修正；分析结果本轮未重新验证。

## 修正决策账本

| 日期 | 现象 | 根因证据 | 最小修正 | 再验证结果 | 是否需人工决定 |
|---|---|---|---|---|---|
| 2026-08-27 | 迭代未覆盖最终 UMAP 优化 | 用户明确要求审查并自主调整降维/聚类参数 | 扩展 README 迭代规则与权限 | Markdown 与规则一致性检查通过；分析未运行 | 否 |
| 2026-08-27 | 文档不支持闭环迭代 | 缺少验证层级、停止条件和记录模板 | 仅重构 README/Report | 文档检查通过；分析未运行 | 否 |

## 新一轮记录模板

复制以下小节到“迭代记录”顶部，并填写实际证据；删除不适用项，不保留占位语句。

```markdown
### YYYY-MM-DD：<本轮目标>

- 请求与权限：<允许读取、执行、修改或覆盖的范围>
- Notebook SHA-256 / Git 状态：<值>
- 环境：<kernel、Python、Scanpy、关键依赖>
- 输入：<路径、身份、维度、校验信息>
- 现有输出：<文件及是否允许覆盖>
- 执行结果：<完成位置、错误、耗时和资源信息>
- 五层验证：<执行、结构、分析、注释、输出>
- 最终注释 UMAP 诊断：<碎裂、混合、样本效应、QC 驱动、稀有群、marker/邻域一致性>
- 参数候选与假设：<本轮只改变的参数家族、候选值、预期改善>
- 候选输出目录：<独立且不覆盖正式结果的路径>
- 与上一基线比较：<改善、退化、稳定性和是否接受>
- 与参考基线差异：<差异及其解释，不以匹配为目标>
- 根因证据：<失败时填写>
- 最小修正：<修改内容；未修改则写明>
- 再验证：<从何处重跑、哪些下游检查通过>
- 结论：<通过、失败或阻塞>
- 未决风险与下一步：<需要用户决定的事项>
```

## 更新规则

- 新状态写在顶部，不静默改写历史结论。
- 失败证据、异常 traceback 和差异摘要应保留可定位路径。
- 参数值以 Notebook 与 `confirmed_parameters.json` 为准；报告只记录本轮关键差异。
- 修改完成后检查 Markdown、路径、事实一致性和 Git diff。
- 若新证据推翻旧结论，新增一条更正记录，说明被更正内容、证据和影响范围。
