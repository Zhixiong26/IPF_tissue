# Scanpy 分析报告 / Analysis report

状态 / Status：CYL/ZCP notebook 的 18-cluster 结果已完成 marker/UMAP 复核并同步至 Python/Slurm 脚本；14、16、17 保留复核层级，待正式运行确认。

CYL 和 ZCP 的本地 10x filtered matrix ZIP 已放入 `Data/Matrix/`：CYL 有 4,264 cells，ZCP 有 4,685 cells，均有 38,606 features。`Notebooks/E_CYL_ZCP_scanpy.ipynb` 以现有 `Notebooks/S12-2N.ipynb` 为模板，但不修改原 notebook，也不沿用 S12 的输入路径、cluster 编号或输出。

notebook 保留 `layers['counts']` 中的原始整数 counts，并在 `obs` 中保留 cohort 与 pre-filter QC。候选 QC 阈值、是否回归 covariates、PCA/邻接参数、Leiden resolution 和 cluster-to-cell-type mapping 均须在 notebook 中审核。保存 cell types/H5AD 的 cell 默认被 `ANALYSIS_CONFIRMED = False` 阻止。

`Scripts/01_run_e_scanpy.py` 现执行分样本 QC、Scrublet、singlet 合并、归一化、batch-aware HVG、PCA、Harmony、Harmony KNN、UMAP、Leiden 0.8、marker 检验和人工注释。所有派生表格、图片及 H5AD 写到用户指定的 `Results/Scanpy/` 子目录。

正式脚本会保存 pre-filter QC、Scrublet per-cell 结果、HVG/PCA 表、cluster/sample 构成、完整 marker 表、全局及专项 marker dotplot、稀有群 QC、细胞注释表、UMAP、参数、软件版本与最终 H5AD。当前 18-cluster 注释只在重跑编号完整匹配 0–17 时自动应用；编号变化时保存诊断结果并要求重新审核。

The production workflow now guards the reviewed 18-cluster set (0–17), exports focused epithelial and rare-cluster evidence, and does not transfer labels when the observed Leiden IDs change.

2026-08-25 校对结果：注释配置校验为 `pass`（18 clusters、15 个合并 cell types、每类 4 markers）；两个 Python 入口通过编译与 `--help` 检查，Slurm wrapper 通过 `bash -n`。新增审计表连接逻辑通过 18-cluster 合成数据 smoke test。尚未执行全量表达矩阵流程。

Audit result on 2026-08-25: annotation configuration passed (18 clusters, 15 merged cell types, four markers per type); both Python entry points passed compilation and CLI checks, the Slurm wrapper passed `bash -n`, and the joined audit-table logic passed an 18-cluster synthetic smoke test. The full expression workflow has not been run.

## 正式脚本与运行入口

- Python：`Scripts/Scanpy/Scripts/01_run_e_scanpy.py`
- 注释配置校验 / Annotation configuration validator：`Scripts/Scanpy/Scripts/02_validate_annotation_config.py`
- Slurm：`Scripts/Scanpy/Scripts/run_e_scanpy.sbatch`
- 项目根目录：`/home/lijia/luozhixiong/IPF_tissue`
- 默认输出：`Results/Scanpy/E_CYL_ZCP`
- 输出限制：`--output-dir` 必须位于项目 `Results/` 下。
- 随机种子：`0`，用于 Scrublet、PCA、Harmony、邻接图、UMAP 和 Leiden。
- 当前检查环境：Python 3.9、Scanpy 1.9.3、harmonypy 0.0.10、Scrublet 0.2.3；正式运行时的实际版本另存为 `software_versions.json`。

### 完整命令行参数表

| 参数 | 类型 | 默认值 | 用途 |
|---|---|---|---|
| `--cyl-zip` | Path | `Data/Matrix/25100718_CYL_E/filtered_feature_bc_matrix.zip` | CYL 10x ZIP 输入 |
| `--zcp-zip` | Path | `Data/Matrix/25100718_ZCP_E/filtered_feature_bc_matrix.zip` | ZCP 10x ZIP 输入 |
| `--output-dir` | Path | `Results/Scanpy/E_CYL_ZCP` | 结果根目录 |
| `--min-genes` | int | 200 | 单细胞最低检测基因数 |
| `--max-genes` | int | 6000 | 单细胞检测基因数上限 |
| `--min-counts` | int | 500 | 单细胞最低总 counts |
| `--max-mt-percent` | float | 5.0 | 单细胞线粒体 counts 比例上限 |
| `--min-cells-per-gene` | int | 3 | cohort 内基因最低检出细胞数 |
| `--doublet-rate-per-1000` | float | 0.004 | Scrublet 动态先验系数 |
| `--target-sum` | float | 10000 | 每细胞归一化总量 |
| `--n-hvg` | int | 2000 | 高变基因数量 |
| `--pca-components` | int | 50 | 请求计算的 PCA 成分数 |
| `--n-pcs` | int | 30 | Harmony 邻接图使用的成分数 |
| `--n-neighbors` | int | 15 | KNN 近邻数 |
| `--leiden-resolution` | float | 0.8 | Leiden 聚类分辨率 |
| `--seed` | int | 0 | 全流程随机种子 |
| `--overwrite` | flag | False | 允许覆盖同名输出文件 |
| `--allow-unreviewed-clusters` | flag | False | cluster 编号变化时允许探索性运行正常返回 |

### Slurm 资源

| 参数 | 设置 |
|---|---:|
| job name | `ipf_scanpy_e` |
| partition | `fat` |
| nodes | 1 |
| CPUs | 16 |
| memory | 64 GB |
| wall time | 1 day |
| stdout | `Scripts/Scanpy/logs/scanpy_e_%j.out` |
| stderr | `Scripts/Scanpy/logs/scanpy_e_%j.err` |
| Python | `/home/lijia/jiangyuanpei/miniforge3/envs/allcools/bin/python` |

## 输入和读取参数

| cohort | 默认输入 |
|---|---|
| CYL | `Data/Matrix/25100718_CYL_E/filtered_feature_bc_matrix.zip` |
| ZCP | `Data/Matrix/25100718_ZCP_E/filtered_feature_bc_matrix.zip` |

- 输入为 ZIP 内的 `filtered_feature_bc_matrix/` 10x 稀疏矩阵。
- `sc.read_10x_mtx()` 使用 `var_names='gene_symbols'` 和 `make_unique=True`。
- cell barcode 分别添加 `CYL_` 或 `ZCP_` 前缀，样本来源写入 `obs['cohort']`。
- 当前 `obs['sample']` 由 `cohort` 明确复制，因此 sample 类别为 CYL、ZCP；Harmony 仍使用 `cohort` 作为 batch key。
- 当前没有经过验证的生物学 `group` 元数据。脚本预留 `obs['group']` 检查，但不会用 sample 伪造 group。
- 原始整数矩阵复制到 `layers['counts']`。
- 线粒体基因：基因名转大写后以 `MT-` 开头。
- 核糖体基因：基因名转大写后以 `RPL` 或 `RPS` 开头。
- `calculate_qc_metrics()` 使用 `qc_vars=['mt', 'ribo']`、`percent_top=None`、`log1p=False`、`inplace=True`。

## QC 与 Scrublet 参数

QC 和 doublet 检测均对 CYL、ZCP 分别执行，只有通过 QC 的 singlet 才会合并。

| 参数 | 默认值 | 判定 |
|---|---:|---|
| `--min-genes` | 200 | `n_genes_by_counts >= 200` |
| `--max-genes` | 6000 | `n_genes_by_counts < 6000` |
| `--min-counts` | 500 | `total_counts >= 500` |
| `--max-mt-percent` | 5.0 | `pct_counts_mt < 5.0` |
| `--min-cells-per-gene` | 3 | 每个 cohort 内 `filter_genes(min_cells=3)` |
| `--doublet-rate-per-1000` | 0.004 | 见下式 |

每个 cohort 的 Scrublet 先验概率为：

```text
expected_doublet_rate = 0.004 × 通过基础 QC 的细胞数 / 1000
```

Scrublet 使用 `layers['counts']`，显式参数如下：

- `random_state=0`
- `n_prin_comps=30`
- `use_approx_neighbors=False`，使用精确近邻。
- `verbose=True`
- 其余 Scrublet 参数沿用所记录软件版本的库默认值。
- 若任一 cohort 在基础 QC 后少于 100 个细胞，脚本停止。
- 自动阈值产生的 `predicted_doublet=True` 细胞从下游排除，但分数和预测结果全部保存。

QC 描述统计保存 `count/mean/std/min/1%/5%/50%/95%/99%/max`。每个 cohort 绘制最高表达的 20 个基因、四项 QC 小提琴图、counts–MT 和 counts–genes 散点图，以及 Scrublet observed/simulated 分数图。

## 合并、归一化和 HVG 参数

- 合并对象：各 cohort 通过 QC 且不是 Scrublet doublet 的细胞。
- `ad.concat()`：`join='outer'`、`merge='same'`、`index_unique=None`。
- 合并后检查 cell ID 不重复，并把 `cohort` 转换为 category。
- `normalize_total(target_sum=10000)`。
- `log1p()`。
- 归一化后的全基因对象保存为 `adata.raw`，供 marker 分析和 dotplot 使用。
- `highly_variable_genes()`：`n_top_genes=2000`、`batch_key='cohort'`、`flavor='seurat'`。
- 下游工作对象仅保留 HVG。
- `REGRESS_COVARIATES=False`：不回归 `total_counts` 或 `pct_counts_mt`。
- `scale(max_value=10)`：中心化、标准化并截断极端标准化值。

## PCA 与 Harmony 参数

### PCA

| 参数 | 设置 |
|---|---:|
| 请求 PC 数 | 50 |
| 实际 PC 数 | `min(50, n_cells - 1, n_HVG - 1)` |
| solver | `arpack` |
| random state | 0 |
| 邻接图使用 PC 数 | 30 |

若实际可计算 PC 少于邻接图请求的 30 个，脚本停止。输出每个 PC 的解释方差比例、方差曲线及 Harmony 前按 cohort 着色的 PCA。

### Harmony

| 参数 | 设置 |
|---|---|
| batch key | `cohort` |
| input basis | `X_pca` |
| output basis | `X_pca_harmony` |
| `nclust` | `min(100, max(2, round(n_cells / 30)))` |
| `sigma` | 长度等于 `nclust`、每项均为 `0.1` 的数组 |
| maximum iterations | 20 |
| random state | 0 |

除上述显式设置外，其余 Harmony 参数使用 harmonypy 运行版本的默认值。

## 邻接图、UMAP 与 Leiden 参数

| 步骤 | 参数 |
|---|---|
| Harmony 前 KNN | `use_rep='X_pca'`, `key_added='neighbors_before_harmony'` |
| Harmony 后 KNN | `use_rep='X_pca_harmony'`，作为 Leiden 正式邻接图 |
| 两套 KNN | `n_neighbors=15`, `n_pcs=30`, `random_state=0` |
| Harmony 前 UMAP | `neighbors_key='neighbors_before_harmony'`, `random_state=0` |
| Harmony 后 UMAP | 默认 Harmony 邻接图，`random_state=0` |
| Harmony 前坐标 | `obsm['X_umap_before_harmony']` |
| Harmony 后坐标 | `obsm['X_umap_after_harmony']`，同时保留为默认 `X_umap` |
| Leiden resolution | 0.8 |
| Leiden seed | 0 |

去批次前后使用完全相同的近邻数、PC 数和 UMAP 随机种子，唯一主要差异是输入表示分别为原始 PCA 与 Harmony PCA。Leiden 仅在 Harmony 后正式邻接图上计算。未显式覆盖的 UMAP、KNN 和 Leiden 参数使用 Scanpy 1.9.3 对应默认值。

脚本保存一张按 sample 着色的 Harmony 前后并排对照图，以及 Harmony 后独立的 cell type、sample 和 group UMAP。当前没有 group 元数据，因此不生成虚假的 group 图，而是写出 `group_metadata_status.json` 说明跳过原因；将来加入经验证的 `obs['group']` 后会自动生成独立 group UMAP。

## Marker 检验和注释验证参数

- `rank_genes_groups(groupby='leiden', method='wilcoxon', use_raw=True)`。
- 导出所有 cluster 的完整 marker 统计量及每群前 50 个 marker 名称。
- marker 排名图显示每群前 20 个基因，`sharey=False`。
- dotplot 在绘图前按最终 `cell_type` 合并 cluster；例如 cluster 5/9 合并为一行 Macrophages，cluster 7/8 合并为一行 Endothelial cells。
- 左侧纵轴和顶部 marker 分组使用完全相同的 cell-type 名称集合。两轴以同一个 `CELL_TYPE_ORDER` 初始化；dendrogram 重排时同步重排两轴，保证顺序始终一致。
- 列为 marker genes，并按 cell type 分组；每组配置 4 个候选 marker，自动删除数据中不存在的基因，因此实际显示 3–4 个，顶部显示竖排分组名和括号。
- `use_raw=True`，不做 `standard_scale`；红色色阶直接表示每群的平均 log-normalized expression，点大小表示表达该基因的细胞比例。
- 右侧 dendrogram 使用展示 marker 的群均值、Pearson correlation 和 complete linkage 计算。
- 样式为 `cmap='Reds'`、`dot_min=0`、`dot_max=0.6`、`smallest_dot=0`、`largest_dot=180`、灰色点边框 `0.5`、`size_exponent=1.5`。

### 正式 dotplot 精简 marker（每类显示 3–4 个）

| 顶部分组 | marker |
|---|---|
| AT2 | SFTPC, ABCA3, NAPSA, LPCAT1 |
| Secretory epithelial | NEDD4L, SFTA3, SCNN1B, GPRC5A |
| Fibroblasts | COL1A1, COL1A2, DCN, COL3A1 |
| Ciliated cells | FOXJ1, DNAH11, PIFO, CFAP46 |
| Secretory / mucous epithelial | BPIFB1, MUC4, WFDC2, TMC5 |
| Macrophages | LST1, C1QA, MRC1, CD163 |
| AT1-like | AGER, CAV1, HOPX, AQP5 |
| Endothelial cells | PECAM1, VWF, KDR, EMCN |
| Basal cells | KRT5, KRT15, TP63, KRT17 |
| Smooth muscle / mural cells | ACTA2, MYH11, CARMN, PDGFRB |
| MT-high AT2-like | MT-CO1, MT-ND1, MT-ND4, MT-CYB |
| T cells | CD3D, CD3E, TRAC, BCL11B |
| Cycling cells | MKI67, TOP2A, RRM2, ANLN |
| Lymphatic endothelial cells | PROX1, FLT4, CCL21, LYVE1 |
| NA | EPCAM, KRT8, COL11A1, CEMIP |

The formal dotplots merge Leiden clusters by reviewed cell type before plotting. The y-axis labels and top marker-group labels are the same cell-type set and are reordered together by the marker-expression dendrogram.

### 注释 marker 集合

| 类型 | marker |
|---|---|
| Pan-epithelial | EPCAM, KRT8, KRT18, KRT19 |
| AT2 | SFTPC, SFTPB, SFTPA1, SFTPA2, ABCA3, LPCAT1 |
| AT1 | AGER, CAV1, CAV2, PDPN, HOPX, EMP2, AQP5 |
| Secretory | SCGB1A1, SCGB3A1, SCGB3A2, BPIFB1, MUC4, WFDC2, SLPI |
| Basal | KRT5, KRT14, KRT15, KRT17, TP63, MIR205HG |
| Ciliated | FOXJ1, PIFO, TPPP3, DNAH11, CFAP46, HYDIN |
| Macrophage | LST1, TYROBP, FCER1G, CD68, C1QA, C1QB, C1QC, MRC1, CD163, PPARG |
| Fibroblast | COL1A1, COL1A2, COL3A1, DCN, LUM, PDGFRA, COL6A3 |
| Pan-endothelial | PECAM1, VWF, KDR, EMCN, ENG, ESAM, RAMP2, PLVAP |
| Capillary endothelial | CA4, RGCC, EMCN, GPIHBP1, BTNL9, EDNRB, EPAS1 |
| Lymphatic endothelial | PROX1, PDPN, LYVE1, FLT4, CCL21, MMRN1, RELN |
| Smooth muscle/mural | ACTA2, TAGLN, MYH11, LMOD1, CNN1, CARMN, PDGFRB, RGS5 |
| T cell | CD3D, CD3E, TRAC, BCL11B, ITK, CD247, IL7R |
| Mast | KIT, CPA3, TPSAB1, TPSB2, HDC, MS4A2, HPGDS |
| Cycling | MKI67, TOP2A, UBE2C, CENPF, BIRC5, RRM2, ANLN, ECT2, DIAPH3 |

### 人工审核的 cluster 注释

| Leiden | cell type | confidence | 主要证据 |
|---:|---|---|---|
| 0 | AT2 | high | SFTPC/SFTPB/ABCA3/LPCAT1 |
| 1 | Secretory epithelial | medium | NEDD4L/SFTPB/SFTA3；位于 epithelial compartment |
| 2 | Fibroblasts | high | COL1A2/COL3A1/COL5A1/COL6A3/PDGFRA |
| 3 | Ciliated cells | high | CFAP/DNAH/HYDIN |
| 4 | Secretory / mucous epithelial | high | BPIFB1/MUC4/ERN2/TMC5 |
| 5 | Macrophages | high | CD163/MRC1/CTSB/FCER1G |
| 6 | AT1-like | medium | CAV1/HOPX/CAV2 上升，但 AGER/PDPN/AQP5 不完整且保留 SFTPB/LPCAT1/ABCA3 |
| 7 | Endothelial cells | high | EPAS1/PECAM1/VWF/BTNL9 |
| 8 | Endothelial cells | high | VWF/PTPRB/PECAM1/EPAS1 |
| 9 | Macrophages | high | PPARG/MRC1/CD163/MSR1 |
| 10 | Basal cells | high | EGFR/KRT15/COL7A1/TP63/KRT5 |
| 11 | Smooth muscle / mural cells | high | MYH11/LMOD1/CARMN/PDGFRB/COL4A1 |
| 12 | MT-high AT2-like | medium | MT genes 与 SFTPC/SFTPA2/SFTPB/ABCA3；状态标签 |
| 13 | T cells | high | PTPRC/CD247/BCL11B/ITK/DOCK2 |
| 14 | Cycling cells | medium | DIAPH3/FANCI/MELK/RRM2/ECT2/ANLN；需泛上皮 panel 确认谱系 |
| 15 | Lymphatic endothelial cells | high | PROX1/FLT4/LYVE1/CCL21/MMRN1/RELN 与 pan-endothelial 信号 |
| 16 | NA | low | 50 cells，94% CYL；保留 unresolved epithelial-like 证据备注 |
| 17 | NA | low | 20 cells，90% ZCP；保留 COL1A2/DCN/COL3A1 stromal-like 证据备注 |

### 逐 cluster 注释验证记录 / Per-cluster annotation validation

下表基于当前 resolution 0.8 的 18-cluster 结果、每群 Top markers、经典 marker 表达比例、QC/样本构成、UMAP 邻接关系及既往粗粒度注释综合审核。细胞数和 CYL/ZCP 比例是本次验证运行的观察值；正式重跑后必须重新核对，不应脱离当前 cluster 集合复用。

| Leiden | Cells；CYL/ZCP | 主要阳性证据 | 排除项、限制及交叉验证 | 审核结论 |
|---:|---|---|---|---|
| 0 | 1261；62.7%/37.3% | SFTPC 79.6%、SFTPB 99.8%、ABCA3 87.5%、LPCAT1 84.1%；Top markers 含 SFTPB/ABCA3 | AT1、basal、secretory 和 immune 程序不占主导；与既往 AT2 一致 | **AT2，high** |
| 1 | 1014；57.0%/43.0% | NEDD4L、SFTA3、GPRC5A/MAGI3，并保留 SFTPB 97.4% | 经典 SCGB1A1/BPIFB1 较弱，且带明显 AT2-like 表达，因此不定为成熟 goblet；与旧 Goblet/Secretory 大类部分对应 | **Secretory epithelial，medium** |
| 2 | 993；39.1%/60.9% | LAMA2、COL6A3、COL1A2、COL3A1、PDGFRA/LUM/DCN extracellular-matrix 程序 | 缺少 epithelial、endothelial 和 immune 主导程序；与既往 Fibroblast 一致 | **Fibroblasts，high** |
| 3 | 982；24.0%/76.0% | CFAP46/54/100、DNAH3/10/11/12、HYDIN、PIFO/TPPP3 纤毛轴丝程序 | 其他 lineage marker 仅为背景/ambient；与既往 Ciliated 一致 | **Ciliated cells，high** |
| 4 | 684；34.1%/65.9% | BPIFB1 68.3%、MUC4、TMC5、ERN2、CHST9；SCGB1A1/MUC5AC 有限但可检出 | 与 cluster 1 同属旧 Goblet/Secretory 大类，但 mucous/airway-secretory 程序更明确 | **Secretory / mucous epithelial，high** |
| 5 | 675；35.0%/65.0% | CD163、MRC1、CTSB、FMNL2、FCER1G/ITGAX 髓系程序 | 与 cluster 9 邻近但 PPARG/MRC1 状态不同；与旧 Macrophage/Myeloid 一致 | **Macrophages，high** |
| 6 | 587；50.1%/49.9% | CAV1 64.9%、HOPX 62.5%、CAV2 36.3%，并有 EMP2/SCEL | AGER 15.0%、PDPN 10.9%、AQP5 2.7% 不完整；仍保留 SFTPB 95.6%、LPCAT1 72.7%、ABCA3 48.7%、SFTPC 33.0%，不支持纯 AT1；旧 AT1 提供方向性支持 | **AT1-like，medium** |
| 7 | 554；61.6%/38.4% | EPAS1、PECAM1、PTPRM、VWF/BTNL9，pan/capillary endothelial 程序 | 与 cluster 8 相邻但 marker state 不同；一级注释与旧 Endothelial 一致 | **Endothelial cells，high** |
| 8 | 367；25.3%/74.7% | VWF、PTPRB、PECAM1、EPAS1、ANO2/SLCO2A1 血管内皮程序 | 与 cluster 7 合并为同一一级 cell type，但保留独立 Leiden state；与旧 Endothelial 一致 | **Endothelial cells，high** |
| 9 | 359；41.5%/58.5% | PPARG、MRC1、MSR1、CD163、SLC11A1、ABCG1 肺巨噬细胞程序 | 与 cluster 5 同谱系但更偏 PPARG/MRC1 resident-like 状态；后续可 subset 重聚类 | **Macrophages，high** |
| 10 | 350；40.0%/60.0% | KRT5 39.7%、KRT15 52.9%、TP63 47.1%、KRT17 25.4%，另有 EGFR/COL7A1/MIR205HG | AT2/secretory marker 不成套；与既往 Basal 一致 | **Basal cells，high** |
| 11 | 197；56.9%/43.1% | MYH11、CARMN、LMOD1、PDGFRB、PRKG1/CALD1 mural-contractile 程序 | 与 endothelial/stromal 区域相邻但缺少 pan-endothelial 主导信号；与旧 Smooth muscle/Pericyte 一致 | **Smooth muscle / mural cells，high** |
| 12 | 193；65.3%/34.7% | Top markers 由 MT-CO/MT-ND/MT-CYB 主导，同时 SFTPC 66.3%、SFTPB 90.7%、ABCA3 61.7%、LPCAT1 70.5% | median pct_counts_mt 约 0.84%，仍低于 5% QC 阈值；这是相对 MT-high 的 AT2 状态，不是独立谱系，也不是自动剔除理由 | **MT-high AT2-like，medium/state** |
| 13 | 189；72.0%/28.0% | PTPRC、SKAP1、CD247、DOCK2、BCL11B、ITK/FYN T-lymphocyte 程序 | 缺少成套 KIT/CPA3/TPSAB1/HDC mast 程序；旧框架未设置 lymphocyte 类，不能以旧 NA 覆盖 | **T cells，high** |
| 14 | 111；78.4%/21.6% | DIAPH3、FANCI/FANCA、MELK、RRM2、ECT2、ANLN/TOP2A 完整 cell-cycle 程序；median genes/counts 为 3763/10089 | 同时有 SFTPB 98.2%、ABCA3 61.3%、LPCAT1 80.2%，提示可能为 cycling epithelial，但一级谱系尚不强制；与旧 Proliferating 一致 | **Cycling cells，medium/state** |
| 15 | 110；47.3%/52.7% | PROX1 29.1%、FLT4 40.0%、LYVE1 21.8%、CCL21 19.1%、MMRN1 59.1%、RELN 49.1%，并有 PECAM1 37.3% | DCN 0%，COL1A1/2/3 阳性比例低，不支持 fibroblast；两个样本均存在，不像单样本伪群；旧框架缺失 lymphatic endothelial 类 | **Lymphatic endothelial cells，high** |
| 16 | 50；94.0%/6.0% | SFTPB 66%、HOPX 58%、CAV1 50%、DNAH11 56%，兼有少量 epithelial/AT1-like 线索 | 无完整 MKI67/TOP2A/UBE2C/CENPF/BIRC5/RRM2/ANLN/ECT2 cycling signature；细胞少且强 CYL 偏倚，marker 混合，按用户指定不强制亚型 | **NA，low**；证据备注为 unresolved epithelial-like |
| 17 | 20；10.0%/90.0% | COL11A1、CEMIP、THSD4、GPC6、BNC2；COL1A2/DCN/COL3A1 分别约 40%/30%/25% | 仅 20 cells、强 ZCP 偏倚，证据支持 stromal 方向但不足以稳定命名独立类型；按用户指定保守处理 | **NA，low**；证据备注为 stromal-like |

逐群表中的“high/medium/low”是人工注释置信度，不是统计显著性。Cluster 5/9、7/8 在一级 cell-type dotplot 中分别合并为 Macrophages、Endothelial cells，但 Leiden 编号、marker 表和每群 QC 仍独立保留；clusters 16/17 在一级图中合并为 NA，但不能据此认为两群具有相同生物学状态。

该映射只适用于已审核的完整 cluster 集合 `0–17`。脚本比较实际 Leiden 编号与该集合：若缺少旧 cluster 或出现新 cluster，未审核编号标为 `Unassigned`，保存 `annotation_mismatch.json` 和全部诊断结果，并默认以非零状态退出。只有人工审核后才可更新映射；`--allow-unreviewed-clusters` 仅允许保存探索性结果，不代表注释已确认。

专项验证固定检查上皮群 `1/4/6/10/12/14/16`，以及稀有群 `15/17`。验证运行中 cluster 15 有 110 cells；其 PROX1、FLT4、LYVE1、CCL21、MMRN1、RELN 阳性比例分别约为 29%、40%、22%、19%、59%、49%，同时 DCN 为 0%，因此以联合证据注释为 lymphatic endothelial，而不是仅凭 UMAP 位置判定。Cluster 14 有 111 cells 并具有完整 cell-cycle marker；在新的泛上皮面板得到正式结果前保留较宽的 `Cycling cells`。Cluster 16/17 因细胞少且样本偏倚明显，按用户指定合并为一级类型 `NA`，但各自的 epithelial-like/stromal-like 线索继续保存在证据列中。

Focused validation covers epithelial clusters `1/4/6/10/12/14/16` and rare clusters `15/17`. Cluster 15 is called lymphatic endothelial from the combined lymphatic and pan-endothelial program, whereas clusters 14, 16, and 17 retain explicit review/low-confidence status.

## 与既往粗粒度注释交叉验证 / Cross-validation against the prior broad annotation

用户提供的同批 CYL/ZCP 肺组织既往结果包含 AT1、AT2、Ciliated、Basal、Goblet/Secretory、Proliferating、Macrophage/Myeloid、Endothelial、Smooth muscle/Pericyte、Fibroblast 和 NA。它作为独立的谱系级辅助证据使用，不作为当前 cell labels 的替代来源。

| 既往类型 | 当前 Leiden | 当前注释 | 判断 |
|---|---|---|---|
| AT2 | 0 | AT2 | 一致 |
| AT1 | 6 | AT1-like | 部分一致；保留过渡标签 |
| Ciliated | 3 | Ciliated cells | 一致 |
| Basal | 10 | Basal cells | 一致 |
| Goblet/Secretory | 1, 4 | 两个 secretory epithelial 群 | 一致并细化 |
| Macrophage/Myeloid | 5, 9 | Macrophages | 一致并保留两个状态 |
| Endothelial | 7, 8 | Endothelial cells | 一致 |
| Smooth muscle/Pericyte | 11 | Smooth muscle / mural cells | 一致 |
| Fibroblast | 2, 17 | Fibroblasts；NA | cluster 17 改为 NA，基质证据仅作备注 |
| Proliferating | 14 | Cycling cells | 一致 |
| NA / rare | 12, 13, 15, 16, 17 | MT-high AT2-like、T cells、Lymphatic endothelial、NA | clusters 16/17 使用 NA；其余群提供额外信息 |

Cluster 6 的扩展 AT1 检查显示：CAV1 约 64.9%、HOPX 62.5%、CAV2 36.3%，但 AGER 15.0%、PDPN 10.9%、AQP5 2.7%；同时 SFTPB 95.6%、LPCAT1 72.7%、ABCA3 48.7%、SFTPC 33.0%。因此旧结果支持其位于 AT1 谱系方向，但当前证据不支持改成纯 `AT1`，继续使用 `AT1-like`。

旧框架没有 T cells 和 lymphatic endothelial 类别，不能用其 NA 标签覆盖当前 cluster 13/15 的明确 marker 证据。Cluster 15 仍保留 `Lymphatic endothelial cells`；cluster 13 仍保留 `T cells`。对应关系保存于 `Supplementary/prior_annotation_crosswalk.tsv`。

The prior result supports the major lineage structure but is intentionally treated as auxiliary evidence. Because only the published/illustrated broad labels—not per-cell barcode-level annotations—are currently available, this is a lineage crosswalk rather than a quantitative cell-by-cell confusion matrix.

## 图片参数

- Matplotlib 非交互后端：`Agg`。
- Scanpy 全局：`dpi=100`、`frameon=False`。
- 保存 PNG：`dpi=220`、`bbox_inches='tight'`、白色背景。
- Harmony 前后 sample 对照：两个独立坐标面板，`size=10`、`alpha=0.75`，图幅 `13 × 5`。
- cell-type UMAP：Harmony 后坐标，`size=10`、`alpha=0.85`、`legend_loc='right margin'`。
- sample UMAP：Harmony 后坐标，`size=10`、`alpha=0.75`、`legend_loc='right margin'`。
- group UMAP：仅在真实 `obs['group']` 存在且至少有一个非缺失值时绘制；参数与 sample UMAP 相同。
- sample 配色：CYL `#0072B2`，ZCP `#E69F00`。
- Leiden UMAP：数字标签 `legend_loc='on data'`。

## 输出文件和 H5AD 结构

结果目录下保存：

- `parameters.json`：Python 脚本本次实际参数，包括 marker 检验方法、审核 cluster 集合、专项 cluster 列表与全部 marker panels。
- `software_versions.json`：运行环境版本。
- `run_summary.json`：输入维度、各阶段细胞数、cluster 数和注释状态。
- `annotation_guard_status.json`：每次运行均保存实际与审核 Leiden 集合的精确比较。
- `annotation_mismatch.json`：仅在 cluster 编号与审核结果不一致时产生。
- `group_metadata_status.json`：记录 group 元数据是否可用以及 group UMAP 是否生成。
- `tables/cell_qc_prefilter.tsv.gz`：所有输入细胞的 QC 与基础过滤标签。
- `tables/qc_distribution_summary.tsv`：分 cohort QC 描述统计。
- `tables/qc_doublet_summary.tsv`：QC、doublet 与 singlet 数量汇总。
- `tables/scrublet_per_cell.tsv.gz`：进入 Scrublet 的每个细胞及其分数、预测。
- `tables/highly_variable_genes.tsv`：HVG 指标。
- `tables/pca_variance_ratio.tsv`：每个 PC 的解释方差。
- `tables/cluster_sample_counts.tsv` 和 `cluster_sample_fractions.tsv`。
- `tables/cluster_qc_summary.tsv`：每群 QC 中位数。
- `tables/rare_cluster_qc_summary.tsv`：cluster 15/17 的细胞数、样本构成、QC 与 Scrublet 分数。
- `tables/leiden_ranked_markers.tsv.gz`：完整 marker 统计。
- `tables/leiden_top50_marker_names.tsv`：每群前 50 marker。
- `tables/annotation_marker_presence.tsv`：经典 marker 是否存在于数据中。
- `tables/dotplot_marker_panel.tsv`：合并 cell type、分组/marker 顺序及 marker 是否存在。
- `figures/annotation/epithelial_focus_marker_dotplot.png`：上皮 compartment 与 cycling 专项验证。
- `figures/annotation/rare_cluster_marker_dotplot.png`：cluster 15/17 的谱系排查。
- `tables/reviewed_cluster_annotations.tsv`：人工注释、置信度和证据。
- `tables/cluster_annotation_validation.tsv`：逐 cluster 注释证据、是否在本次运行出现、QC 和 sample 构成的合并审计表。
- `tables/cell_type_cluster_membership.tsv`：合并后 cell type 与成员 Leiden clusters 的显式对应关系。
- `tables/cell_type_sample_counts.tsv` 和 `cell_type_sample_fractions.tsv`。
- `tables/cell_metadata_qc_annotation.tsv.gz`：最终每细胞 QC、Scrublet、cohort、sample、可选 group、Leiden、cell type，以及 Harmony 前后两套 UMAP 坐标。
- `figures/qc/`、`figures/doublet/`、`figures/dimensionality/`、`figures/clustering/`、`figures/annotation/`：对应分析阶段的 PNG。
- `objects/rna_e_cyl_zcp_annotated.h5ad`：cluster 与审核集合一致时的正式对象。
- `objects/rna_e_cyl_zcp_unreviewed_clusters.h5ad`：编号变化时的待审核对象。

最终 H5AD 的主矩阵 `X` 为缩放后的 2000 HVG，`raw` 为合并 singlet 的全基因 log-normalized 表达，`layers['counts']` 在 HVG 工作对象中保留对应 HVG 的原始 counts。`obsm` 保存原始 PCA、Harmony PCA、Harmony 前 UMAP、Harmony 后 UMAP、默认正式 UMAP 及兼容别名 `X_umap_harmony`；`obs` 保存 cohort、sample、可选 group、QC、Scrublet、Leiden 和 cell type；`uns` 保存运行参数与注释证据。

`02_validate_annotation_config.py` 仅校验生产脚本中注释配置的结构一致性，不验证表达数据中的生物学身份。正式结果仍须通过 `annotation_guard_status.json`、逐群 marker、QC/sample 构成和专项 dotplot 复核。若实际 cluster 集合不完全等于 0–17，所有细胞统一标为 `Unassigned`；即使某个数字编号仍存在，也不会沿用旧标签。逐群审计表会同时保留缺失的审核群和新出现的未审核群。

`02_validate_annotation_config.py` checks structural consistency of the annotation configuration only; it does not biologically validate expression-derived identities. A formal run still requires review of the annotation guard, cluster markers, QC/sample composition, and focused dotplots. If the observed cluster set is not exactly 0–17, every cell remains `Unassigned`; numeric IDs that happen to persist do not inherit old labels. The per-cluster audit retains both missing reviewed clusters and newly observed clusters.

## 命令行参数和覆盖策略

所有 CLI 默认值已在上述章节记录。附加开关：

- `--overwrite`：允许覆盖已存在的同名输出；不会自动删除目录内其他旧文件，因此正式分析仍建议使用新的空目录。
- `--allow-unreviewed-clusters`：cluster 编号变化时允许脚本正常返回，但输出仍标记为 `annotation_requires_review`。
- Slurm wrapper 默认拒绝非空输出目录，也拒绝 `Results/` 以外的路径。

示例提交：

```bash
cd /home/lijia/luozhixiong/IPF_tissue
sbatch Scripts/Scanpy/Scripts/run_e_scanpy.sbatch \
  /home/lijia/luozhixiong/IPF_tissue/Results/Scanpy/E_CYL_ZCP_YYYYMMDD
```
