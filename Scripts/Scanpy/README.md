# Scanpy

本阶段保存针对明确 AnnData 输入的版本化预处理、质控、邻接图、聚类、注释和可视化脚本。

This stage contains versioned preprocessing, QC, neighbourhood, clustering, annotation, and visualization scripts operating on explicit AnnData inputs.

原始输入从项目 `Data/Matrix/` 读取，派生结果写入 `Results/Scanpy/`，运行日志写入 `Scripts/Scanpy/logs/`。正式脚本记录 H5AD 结构、细胞 QC、Scrublet、Harmony 表示、聚类、注释证据、软件环境和随机种子。

Read immutable inputs below `Data/Matrix/`, write derived outputs below `Results/Scanpy/`, and write runtime logs below `Scripts/Scanpy/logs/`.

## Directory layout

- `Notebooks/`: interactive analysis notebooks and the unchanged reference notebook.
- `Scripts/01_run_e_scanpy.py`: notebook-confirmed production workflow.
- `Scripts/02_validate_annotation_config.py`: fast structural check for the reviewed cluster mapping, merged cell types, and dotplot marker panel.
- `Scripts/run_e_scanpy.sbatch`: Slurm production wrapper.
- `logs/`: Slurm standard output and error logs.
- `Report.md`: current analysis status and decisions.

## Transcriptome E: CYL + ZCP

The canonical analysis is notebook-first. Start with `Notebooks/E_CYL_ZCP_scanpy.ipynb`, adapted from `Notebooks/S12-2N.ipynb`, and confirm QC thresholds, representation, clustering, markers and annotation interactively. The original S12 notebook remains unchanged.

`Scripts/01_run_e_scanpy.py` is the production conversion of the confirmed notebook workflow. It performs per-cohort QC and Scrublet filtering, merges singlets, selects batch-aware HVGs, runs PCA and Harmony, constructs matched before/after-Harmony neighbour graphs and UMAPs, computes Leiden markers on the Harmony graph, validates the reviewed 18-cluster mapping, and exports global plus focused epithelial/rare-cluster evidence. Marker dotplots merge clusters by reviewed cell type, display three to four available genes per cell type, and use the same dendrogram-synchronized cell-type order on the left axis and top marker groups. Final UMAPs are separate for cell type and sample; a third group UMAP is generated only after verified `group` metadata is supplied.

The notebook keeps threshold candidates in one cell and prevents final output saving until `ANALYSIS_CONFIRMED = True`. The workflow is unsupervised; methylation-derived manual cell types and S12 cluster-number mappings are not transferred to CYL/ZCP RNA barcodes automatically.

正式运行前先执行轻量注释配置检查。它不读取表达矩阵，也不替代生物学复核；它检查 0–17 映射是否连续、置信度是否合法、合并后的 cell types 是否与点图分组一致，以及每类是否保留 3–4 个不重复 marker。

Before a formal run, execute the lightweight annotation-configuration check. It does not read the expression matrices or replace biological review; it verifies cluster IDs, confidence values, merged cell-type groups, and the three-to-four-marker dotplot contract.

```bash
cd /home/lijia/luozhixiong/IPF_tissue
/home/lijia/jiangyuanpei/miniforge3/envs/allcools/bin/python \
  Scripts/Scanpy/Scripts/02_validate_annotation_config.py
```

A successful report has `"status": "pass"` and explicitly lists merged memberships such as Macrophages = 5,9; Endothelial cells = 7,8; and NA = 16,17.

Run on Slurm with an empty result directory:

```bash
mkdir -p /home/lijia/luozhixiong/IPF_tissue/Scripts/Scanpy/logs
cd /home/lijia/luozhixiong/IPF_tissue
sbatch Scripts/Scanpy/Scripts/run_e_scanpy.sbatch \
  /home/lijia/luozhixiong/IPF_tissue/Results/Scanpy/E_CYL_ZCP_YYYYMMDD
```

The result directory contains `tables/`, `objects/`, and figure subdirectories for QC, Scrublet, dimensionality reduction, clustering, and annotation. Annotation auditing includes `annotation_guard_status.json`, `tables/cluster_annotation_validation.tsv`, `tables/cell_type_cluster_membership.tsv`, and `tables/dotplot_marker_panel.tsv`. If a fresh run does not reproduce the reviewed cluster IDs 0–17, every cell remains `Unassigned`, diagnostic outputs are retained, and the job exits non-zero; review all clusters before updating the mapping.
