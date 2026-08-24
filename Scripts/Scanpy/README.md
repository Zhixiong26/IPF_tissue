# Scanpy

本阶段保存针对明确 AnnData 输入的版本化预处理、质控、邻接图、聚类、注释和可视化脚本。

This stage contains versioned preprocessing, QC, neighbourhood, clustering, annotation, and visualization scripts operating on explicit AnnData inputs.

输出写入项目 `Data/Scanpy/`，运行日志写入 `Scripts/Scanpy/logs/`。实现流程时记录 H5AD 结构、观测元数据、layers、邻接图使用的表示、软件环境、随机种子和脚本顺序。

Write outputs below the project `Data/Scanpy/` directory and runtime logs below `Scripts/Scanpy/logs/`. Document the H5AD schema, observation metadata, layers, neighbour representation, environment, random seed, and script order.

## Transcriptome E: CYL + ZCP

`01_run_e_scanpy.py` loads the local 10x filtered matrices from `Data/Matrix/`, prefixes barcodes with cohort, records raw integer counts in `layers['counts']`, calculates basic RNA QC, and performs normalization, HVG selection, PCA, neighbours, Leiden clustering, and UMAP.

The default QC thresholds are intentionally explicit and adjustable: at least 200 detected genes, at least 500 counts, at most 25% mitochondrial counts, and genes detected in at least three cells. The workflow is unsupervised; methylation-derived manual cell types are not transferred to RNA barcodes automatically.

```bash
cd /home/lijia/luozhixiong/IPF_tissue
OUTPUT_DIR=$PWD/Data/Scanpy/E_CYL_ZCP_20260824
sbatch Scripts/Scanpy/run_e_scanpy.sbatch "$OUTPUT_DIR"
```

The job writes `rna_e_scanpy.h5ad`, `cell_qc.tsv`, `run_summary.json`, and QC/UMAP figures under `OUTPUT_DIR`.
