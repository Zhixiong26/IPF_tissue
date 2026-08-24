# Scanpy

本阶段保存针对明确 AnnData 输入的版本化预处理、质控、邻接图、聚类、注释和可视化脚本。

This stage contains versioned preprocessing, QC, neighbourhood, clustering, annotation, and visualization scripts operating on explicit AnnData inputs.

输出写入项目 `Data/Scanpy/`，运行日志写入 `Scripts/Scanpy/logs/`。实现流程时记录 H5AD 结构、观测元数据、layers、邻接图使用的表示、软件环境、随机种子和脚本顺序。

Write outputs below the project `Data/Scanpy/` directory and runtime logs below `Scripts/Scanpy/logs/`. Document the H5AD schema, observation metadata, layers, neighbour representation, environment, random seed, and script order.

## Transcriptome E: CYL + ZCP

The canonical analysis is notebook-first. Start with `E_CYL_ZCP_scanpy.ipynb`, adapted from `S12-2N.ipynb`, and confirm QC thresholds, representation, clustering, markers and annotation interactively. The original S12 notebook remains unchanged.

`01_run_e_scanpy.py` and `run_e_scanpy.sbatch` are provisional drafts created before notebook validation. Do not submit them as the final analysis. Update them only after the notebook parameters and outputs are confirmed.

The notebook keeps threshold candidates in one cell and prevents final output saving until `ANALYSIS_CONFIRMED = True`. The workflow is unsupervised; methylation-derived manual cell types and S12 cluster-number mappings are not transferred to CYL/ZCP RNA barcodes automatically.

After notebook validation, convert the confirmed cells into the versioned Python entry point, update the Slurm wrapper, run smoke validation, and then submit the formal job.
