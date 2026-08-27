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
run_methscan.sbatch
run_methscan_common.sbatch
run_methscan_qc_stage.sbatch
run_methscan_branch.sbatch
run_methscan_summary.sbatch
submit_methscan_pipeline.sh
```

一次提交会依次运行 selection、ALLC→COV、prepare、filter、smooth，并在之后并行完成三个 threshold 分支：`0.01`、`0.02`、`0.05`。ALLC→COV 转换使用 `ProcessPoolExecutor`，worker 默认自动等于 `SLURM_CPUS_PER_TASK`，也可用 `IPF_METHSCAN_COV_WORKERS` 覆盖。

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

该命令提交公共作业、prepare、filter、smooth、3 个 threshold 分支和汇总作业；prepare/filter/smooth 各申请 4 CPU/16G，并通过 `afterok` 自动串联。

默认输入为 `Data/ALLCools` 和 `Results/Scanpy/E_CYL_ZCP_notebook/cell_id_cell_type.tsv`。可通过 `00_methscan_config.sh` 中的环境变量覆盖路径或参数。

## 执行阶段

1. `01_select_scanpy_cells.py`：发现 CYL/ZCP ALLC，与 RNA `cell_id` 匹配；排除 `NA`、空值和未注释细胞。
2. `02_convert_allc_to_cov.py`：将入选 ALLC 转为 CpG-only Bismark `.cov.gz`，写入 conversion QC。
3. `03_prepare_methscan.py`：调用 `methscan prepare --input-format bismark`。
4. `methscan filter`、`smooth`：使用默认 QC 和平滑参数。
5. `methscan scan → matrix → 04_vmr_scanpy.py`：对 0.01/0.02/0.05 三个 threshold 分支分别运行。
6. `05_summarize_run.py`：验证三个分支的 cell ID、矩阵和 Scanpy 输出并生成 summary。

## 输出结构

```text
<run>/
├── 00_scanpy_selected/
├── 01_cov/cov/*.cov.gz
├── 02_prepared/
├── 03_filtered/
├── 04_scan/var_0.01|var_0.02|var_0.05/VMRs.bed
├── 05_matrix/var_0.01|var_0.02|var_0.05/
├── 06_scanpy/var_0.01|var_0.02|var_0.05/
├── stage_status.tsv
├── <stage>.resources.json
├── run_summary.json
└── run_summary.tsv
```

完成判据：`run_summary.json` 的 `status=complete`，三个分支的 h5ad、UMAP、matrix 均存在且 cell ID 校验通过。取消或失败的目录保留用于诊断，不覆盖后重跑。

## 独立 smoke test

Smoke test 不调用正式 `run_methscan.sbatch`，只在 `/tmp` 取 CYL/ZCP 各 10 个细胞，依次运行全部阶段和三个 threshold；因此不会写入 `Results/Methscan`。推荐用于改脚本后的回归检查。检查重点包括 ALLC→cov 转换、prepare 的 cell ID 校验、三个分支的矩阵/Scanpy 输出和最终 summary。
