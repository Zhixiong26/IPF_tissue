# MethSCAn workflow（整理版）

本目录按执行顺序编号，解决旧版两个 `02` 文件名冲突：

```text
00_methscan_config.sh
01_select_scanpy_cells.py
02_convert_allc_to_cov.py
03_prepare_methscan.py
04_vmr_scanpy.py
05_summarize_run.py
06_run_with_resources.py
07_methdiff_celltype.py
08_hypo_dmr_heatmaps.py
run_methscan.sbatch
run_methscan_common.sbatch
run_methscan_qc_stage.sbatch
run_methscan_branch.sbatch
run_methscan_methdiff.sbatch
run_methscan_hypo_heatmaps.sbatch
run_methscan_pooled_methdiff.sbatch
run_methscan_summary.sbatch
submit_methscan_pipeline.sh
submit_methscan_pooled_dmr.sh
```

一次提交会依次运行 selection、ALLC→COV、prepare、filter、smooth；smooth 后在 CYL、ZCP 各自样本内执行所有保留 cell type 的两两 `methscan diff`，再完成 raw-p fallback、Top200 hypo-DMR、单细胞矩阵和三套热图。同时并行完成三个 VMR threshold 分支：`0.01`、`0.02`、`0.05`。ALLC→COV 转换使用 `ProcessPoolExecutor`，worker 默认自动等于 `SLURM_CPUS_PER_TASK`，也可用 `IPF_METHSCAN_COV_WORKERS` 覆盖。

hypo-DMR 下游口径借鉴 `scLC_ICI_PBMC` 提交 `0c5c46c101fd669edd8b5c6e4e7dccb17d1f4bc3` 所记录的流程。该项目中后补的“08 raw-p fallback”在逻辑上位于 pairwise diff 与原“05 Top200 / 06 matrix”之间；本项目按真实依赖顺序统一写入 `08_rawp_fallback`、`09_top200_hypo_DMRs`、`10_single_cell_DMR_matrix`、`11_hypo_DMR_heatmaps`。

推荐使用一键提交器，将公共阶段、三个 threshold 分支和最终汇总组成 Slurm 依赖 DAG；三个分支在公共阶段完成后并行运行。原 `run_methscan.sbatch` 仍保留为单作业兼容入口。

## 运行

ALLC 必须已解压到 `Data/ALLCools`；脚本不会解压、复制或修改源文件。输出目录必须是不存在的新目录：

```bash
sbatch --partition=cpu --nodelist=cu03 --cpus-per-task=55 --mem=250G \
  /home/lijia/luozhixiong/IPF_tissue/Scripts/Methscan/run_methscan.sbatch \
  /home/lijia/luozhixiong/IPF_tissue/Results/Methscan/<run_name>
```

并行分支提交：

```bash
/home/lijia/luozhixiong/IPF_tissue/Scripts/Methscan/submit_methscan_pipeline.sh \
  /home/lijia/luozhixiong/IPF_tissue/Results/Methscan/<run_name>
```

该命令提交公共作业、prepare、filter、smooth、pairwise meth-diff、hypo-DMR/矩阵/绘图、3 个 threshold 分支和汇总作业；各阶段通过 `afterok` 自动串联。

对已经完成公共过滤和平滑的运行目录，可另外提交 CYL+ZCP 合并细胞类型 DMR 分支：

```bash
/home/lijia/luozhixiong/IPF_tissue/Scripts/Methscan/submit_methscan_pooled_dmr.sh \
  /home/lijia/luozhixiong/IPF_tissue/Results/Methscan/<run_name>
```

该分支将 CYL、ZCP 中同名 RNA cell type 的细胞合并为 `CYL_ZCP`，然后进行 14 种细胞类型的两两 `methscan diff`。它只计算 DMR，结果独立写入 `<run>/pooled/07_methdiff/`，不会覆盖原有分样本 DMR，也不会运行 Top200、单细胞矩阵或热图。

正式 pooled 作业 `309188` 已完成：6,264 个过滤后细胞、14 种细胞类型、91/91 个比较成功，0 个失败、0 个 fallback，共输出 6,437,767 行 DMR。该 pooled DMR 后续由 `Scripts/Methylvi/vmr_dmr/` 使用 `raw p < 0.01`、`|methdiff| >= 0.25` 的全量 unique-DMR 口径构建 VMR+DMR MethylVI 输入；这不会改变本目录原有的分样本 Top200 热图口径。

默认输入为 `Data/ALLCools` 和 `Results/Scanpy/E_CYL_ZCP_notebook/cell_id_cell_type.tsv`。可通过 `00_methscan_config.sh` 中的环境变量覆盖路径或参数。

## 执行阶段

1. `01_select_scanpy_cells.py`：发现 CYL/ZCP ALLC，与 RNA `cell_id` 匹配；排除 `NA`、空值和未注释细胞。
2. `02_convert_allc_to_cov.py`：将入选 ALLC 转为 CpG-only Bismark `.cov.gz`，写入 conversion QC。
3. `03_prepare_methscan.py`：调用 `methscan prepare --input-format bismark`。
4. `methscan filter`、`smooth`：使用默认 QC 和平滑参数。
5. `07_methdiff_celltype.py`：先以无数据复制的 hard-link view 限定 chr1–22/X/Y；默认在每个样本内进行 RNA cell type 两两 `methscan diff`，传入 `--pool-samples-as CYL_ZCP` 时则跨样本合并同名细胞类型。默认 `min-cells=6`、bandwidth=2,000、stepsize=1,000、threshold=0.02。正式包装器启用 `--resume`，仅复用身份字段一致、状态完成且 DMR 格式有效的 comparison。
6. `08_hypo_dmr_heatmaps.py`：仅对已知 `calc_fdr` 除零错误执行隔离的 raw-p fallback；按第 10 列确定 hypo 细胞类型，保留 raw p `<0.01`、组间 mean ratio 差值 `≥0.25` 的 DMR，并按差值降序为每种细胞取 Top200 后合并重叠区间。
7. 同一脚本从 `01_cov/cov/*.cov.gz` 计算单细胞×合并 DMR 矩阵。每个值是 DMR 内 unique CpG ratio 的等权算术平均，缺覆盖为 `NA`。
8. 每个样本只保留至少分配到一个自身 hypo-DMR 的 cell type 作为热图行，生成 `mean_ratio.png`、DMR-wise `zscore.png` 和 `zscore_colorbar_compressed.png`。压缩 colorbar 图保留原 z-score 数值，只将颜色范围设为 `[-1,1]`；普通 z-score 图显示范围默认 `[-3,3]`。
9. `methscan scan → matrix → 04_vmr_scanpy.py`：对 0.01/0.02/0.05 三个 VMR threshold 分支分别运行。
10. `05_summarize_run.py`：验证 pairwise meth-diff、hypo-DMR 图和三个 VMR 分支后生成 summary。

## 输出结构

```text
<run>/
├── 00_scanpy_selected/
├── 01_cov/cov/*.cov.gz
├── 02_prepared/
├── 03_filtered/
├── 07_methdiff/samples/<sample>/comparisons/<cell_type_A>_vs_<cell_type_B>/DMRs.bed
├── 07_methdiff/pairwise_summary.tsv
├── 07_methdiff/pairwise_summary.json
├── 08_rawp_fallback/fallback_summary.json
├── 09_top200_hypo_DMRs/samples/<sample>/
├── 10_single_cell_DMR_matrix/samples/<sample>/single_cell_DMR_mean_unique_CpG_ratio.tsv.gz
├── 11_hypo_DMR_heatmaps/samples/<sample>/mean_ratio.png
├── 11_hypo_DMR_heatmaps/samples/<sample>/zscore.png
├── 11_hypo_DMR_heatmaps/samples/<sample>/zscore_colorbar_compressed.png
├── pooled/07_methdiff/samples/CYL_ZCP/comparisons/
├── pooled/07_methdiff/pairwise_summary.tsv
├── pooled/07_methdiff/pairwise_summary.json
├── 04_scan/var_0.01|var_0.02|var_0.05/VMRs.bed
├── 05_matrix/var_0.01|var_0.02|var_0.05/
├── 06_scanpy/var_0.01|var_0.02|var_0.05/
├── stage_status.tsv
├── <stage>.resources.json
├── run_summary.json
└── run_summary.tsv
```

主流程完成判据：`run_summary.json` 的 `status=complete`，分样本 pairwise/fallback/Top200/matrix/heatmap summary 完成，三个 VMR 分支的 h5ad、UMAP、matrix 均存在且 cell ID 校验通过。pooled 分支独立以 `pooled/methdiff.COMPLETE` 和 `pooled/07_methdiff/pairwise_summary.json` 的 `status=complete` 为准。取消或失败的目录保留用于诊断，不覆盖后重跑。
