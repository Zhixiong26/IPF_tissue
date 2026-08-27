# CYL/ZCP Scanpy Notebook 迭代报告

本文件是可持续更新的运行账本。最新状态放在最前；历史记录保留结论和证据，不复制 Notebook 中的完整实现。稳定规则见 `README.md`。

## 最新状态

- 状态：迭代范围已扩展到最终注释 UMAP 审查及 PCA/Harmony/邻接图/UMAP/Leiden 参数自主优化；本轮未执行分析。
- Canonical Notebook：`Notebooks/E_CYL_ZCP_scanpy.ipynb`
- 本轮 Notebook 行为：未修改、未执行。
- 正式输出目录：`Results/Scanpy/E_CYL_ZCP_notebook/`
- 已知保存开关：Notebook 当前文本中 `ANALYSIS_CONFIRMED = True`、`OVERWRITE_DATA_OUTPUTS = False`；每次运行前必须重新读取确认。
- 已授权范围：可修改并重跑 PCA、Harmony、邻接图、UMAP、Leiden 参数；候选运行不得覆盖正式结果。
- 当前未决事项：下一次实际重放前需要记录 kernel、依赖、输入状态、现有输出清单和最终注释 UMAP 基线。

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
| 执行完整性 | 未验证 | 本轮未运行 Notebook | 使用兼容 kernel 干净重放 |
| 结构完整性 | 未验证 | 仅有历史记录 | 检查 AnnData 字段、维度和 cell ID |
| 分析合理性 | 未验证 | 历史基线为 8,696 cells、18 clusters | 复核 QC、Harmony、Leiden、marker 和最终注释 UMAP |
| 注释完整性 | 未验证 | 历史映射不可自动转移 | 对本轮每个 cluster 重新核验并检查 UMAP/marker/邻域一致性 |
| 输出一致性 | 未验证 | 已知输出接口 | 比较 H5AD、TSV、JSON 和 manifest |

状态只使用：`未验证`、`通过`、`失败`、`阻塞`。没有证据时不得写为“通过”。

## 迭代记录

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
