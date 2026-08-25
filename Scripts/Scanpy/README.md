# Scanpy

本阶段保存针对明确 AnnData 输入的版本化预处理、质控、邻接图、聚类、注释和可视化脚本。

This stage contains versioned preprocessing, QC, neighbourhood, clustering, annotation, and visualization scripts operating on explicit AnnData inputs.

原始输入从项目 `Data/Matrix/` 读取，派生结果写入 `Results/Scanpy/`，运行日志写入 `Scripts/Scanpy/logs/`。正式脚本记录 H5AD 结构、细胞 QC、Scrublet、Harmony 表示、聚类、注释证据、软件环境和随机种子。

Read immutable inputs below `Data/Matrix/`, write derived outputs below `Results/Scanpy/`, and write runtime logs below `Scripts/Scanpy/logs/`.

## Directory layout

- `Notebooks/`: canonical analysis notebooks; `E_CYL_ZCP_scanpy.ipynb` is the preferred execution entry point.
- `Scripts/01_run_e_scanpy.py`: secondary batch conversion; do not use as the default Scanpy entry point.
- `Scripts/02_validate_annotation_config.py`: fast structural check for the reviewed cluster mapping, merged cell types, and dotplot marker panel.
- `Scripts/03_finalize_celltype_annotation.py`: guarded post-run annotation for the 2026-08-25 17-cluster result; it locally reclusters mixed Leiden cluster 9 and writes a non-overwriting reviewed derivative.
- `Scripts/run_e_scanpy.sbatch`: Slurm production wrapper.
- `logs/`: Slurm standard output and error logs.
- `Report.md`: current analysis status and decisions.

## Transcriptome E: CYL + ZCP

The canonical analysis is notebook-first. **When asked to run Scanpy, default to a clean top-to-bottom execution of `Notebooks/E_CYL_ZCP_scanpy.ipynb`, not `Scripts/`.** Save the executed copy and replay summary under a new named `Results/Scanpy/` root, then review QC, clustering, markers and annotation before opening the final save gate. The original S12 notebook remains unchanged.

`Scripts/01_run_e_scanpy.py` is retained as a secondary batch conversion. It may be used when the user explicitly requests script/Slurm execution, for helper validation/export, or after it has demonstrated per-cell equivalence to the notebook in the intended execution environment. Matching parameters and aggregate cell counts alone are insufficient: job `307511` produced 17 clusters while a clean notebook replay reproduced the reviewed 18 clusters.

The notebook keeps threshold candidates in one cell and prevents final output saving until `ANALYSIS_CONFIRMED = True`. The workflow is unsupervised; methylation-derived manual cell types and S12 cluster-number mappings are not transferred to CYL/ZCP RNA barcodes automatically.

The confirmed notebook synchronizes every displayed analysis figure to `Results/Scanpy/E_CYL_ZCP_notebook/figures/` as a 300-dpi PNG and writes `Results/Scanpy/E_CYL_ZCP_notebook/figure_manifest.json`. A figure-only replay preserves the existing H5AD, per-cell annotation table, and confirmed-parameter JSON by default; it first requires the newly reproduced per-cell annotation to match the saved table exactly. Set `OVERWRITE_DATA_OUTPUTS = True` only after explicit review when the formal data products genuinely need replacement.

确认后的 notebook 会把所有在 Jupyter 中展示的分析图同步保存为 300-dpi PNG，位置为 `Results/Scanpy/E_CYL_ZCP_notebook/figures/`，并生成 `figure_manifest.json`。仅补画图片时，默认不覆盖已有 H5AD、逐细胞注释表和参数 JSON；notebook 会先要求本次逐细胞注释与已保存结果完全一致。只有正式数据确需替换且经过明确复核后，才可把 `OVERWRITE_DATA_OUTPUTS` 改为 `True`。

在显式选择 Python/Slurm 辅助入口前，可执行轻量注释配置检查。它不读取表达矩阵，也不替代 notebook 生物学复核；它检查 0–17 映射是否连续、置信度是否合法、合并后的 cell types 是否与点图分组一致，以及每类是否保留 3–4 个不重复 marker。

Before an explicitly selected Python/Slurm helper run, the lightweight annotation-configuration check may be executed. It does not read the expression matrices or replace notebook review; it verifies cluster IDs, confidence values, merged cell-type groups, and the three-to-four-marker dotplot contract.

```bash
cd /home/lijia/luozhixiong/IPF_tissue
/home/lijia/jiangyuanpei/miniforge3/envs/allcools/bin/python \
  Scripts/Scanpy/Scripts/02_validate_annotation_config.py
```

A successful report has `"status": "pass"` and explicitly lists merged memberships such as Macrophages = 5,9; Endothelial cells = 7,8; and NA = 16,17.

Only when script/Slurm execution is explicitly selected, run with an empty result directory:

```bash
mkdir -p /home/lijia/luozhixiong/IPF_tissue/Scripts/Scanpy/logs
cd /home/lijia/luozhixiong/IPF_tissue
sbatch Scripts/Scanpy/Scripts/run_e_scanpy.sbatch \
  /home/lijia/luozhixiong/IPF_tissue/Results/Scanpy/E_CYL_ZCP_YYYYMMDD
```

The result directory contains `tables/`, `objects/`, and figure subdirectories for QC, Scrublet, dimensionality reduction, clustering, and annotation. Annotation auditing includes `annotation_guard_status.json`, `tables/cluster_annotation_validation.tsv`, `tables/cell_type_cluster_membership.tsv`, and `tables/dotplot_marker_panel.tsv`. If a fresh run does not reproduce the reviewed cluster IDs 0–17, every cell remains `Unassigned`, diagnostic outputs are retained, and the job exits non-zero; review all clusters before updating the mapping.

The formal 2026-08-25 run produced clusters 0–16 and correctly stopped at this guard. After marker review, run-specific annotation was finalized without rerunning global Scanpy:

```bash
cd /home/lijia/luozhixiong/IPF_tissue
/home/lijia/jiangyuanpei/miniforge3/envs/allcools/bin/python \
  Scripts/Scanpy/Scripts/03_finalize_celltype_annotation.py
```

The primary per-cell output is `Results/Scanpy/E_CYL_ZCP_20260825/reviewed_annotation_20260825/tables/celltype_annotations.tsv.gz`. The unresolved 70 cells deliberately use the literal label `NA`; when loading with Pandas use `pd.read_csv(path, sep="\t", keep_default_na=False)` so that this label is not converted to a missing value.

2026-08-25 的正式运行得到 0–16，因此按预期触发注释保护。marker 复核后，收尾脚本仅对混合的 Leiden cluster 9 做局部重聚类，并生成本次运行专属的逐细胞注释和 annotated H5AD。70 个未定类型细胞使用字符串 `NA`；Pandas 读取时必须设置 `keep_default_na=False`。

若 reviewed H5AD 已经存在而只需补画最终 UMAP，可运行：

```bash
/home/lijia/jiangyuanpei/miniforge3/envs/allcools/bin/python \
  Scripts/Scanpy/Scripts/03_finalize_celltype_annotation.py \
  --plot-only \
  --publish-figure Results/Scanpy/E_CYL_ZCP_20260825/figures/annotation/umap_cell_type_reviewed.png
```

`--plot-only` 仍会读取原运行的 annotation guard。若原簇集合发生变化，程序会警示该图为 marker 审核后的派生图，并把警示、源 guard 状态、细胞数和输出路径写入 `reviewed_annotation_20260825/umap_cell_type_plot_status.json`。它不会把原始运行重新描述为自动注释成功。

When a reviewed H5AD already exists, `--plot-only` generates the final Harmony-coordinate cell-type UMAP without rerunning annotation. A source guard mismatch is emitted as a warning and preserved in the machine-readable plot-status record; the resulting figure is explicitly a reviewed post-run derivative.
