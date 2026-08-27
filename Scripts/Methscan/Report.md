# MethSCAn workflow report（整理版）

## 脚本编号

旧版 `02_prepare_methscan.py` 与 `02_convert_allc_to_cov.py` 重号，已按实际依赖重排为：

```text
01_select_scanpy_cells.py
02_convert_allc_to_cov.py
03_prepare_methscan.py
04_vmr_scanpy.py
05_summarize_run.py
06_run_with_resources.py
```

主入口 `run_methscan.sbatch` 已同步更新引用，不再调用不存在的 `01_prepare_allc_inputs.py`，也不再使用旧的 `link_path`/`00_manifest` intake 阶段。

新增 `run_methscan_common.sbatch`、`run_methscan_qc_stage.sbatch`、`run_methscan_branch.sbatch`、`run_methscan_summary.sbatch` 和 `submit_methscan_pipeline.sh`。后续正式运行只需执行一次 `submit_methscan_pipeline.sh`：它自动提交公共阶段、三个独立 QC 阶段、三个 threshold 分支和 summary，并设置 `afterok` 前后依赖。prepare/filter/smooth 各申请 4 CPU/16G。

## 当前流程

```text
ALLC → Scanpy whitelist → ALLC→cov → prepare → filter → smooth
→ (scan → matrix → Scanpy) × {0.01, 0.02, 0.05} → summary
```

三套分支共享前面的 prepared/filtered 数据，分别写入 `04_scan`、`05_matrix` 和 `06_scanpy` 下的 `var_0.01`、`var_0.02`、`var_0.05`。

ALLC→cov 阶段按细胞使用多进程转换；`run_methscan.sbatch` 将 worker 数设置为 Slurm 分配的 CPU 数，也可通过 `IPF_METHSCAN_COV_WORKERS` 手动覆盖。当前正式任务在该修改前已使用 16 workers 完成转换；后续新任务会自动跟随申请的 CPU 数。

## 关键筛选规则

- canonical cell ID：`<sample_id>_<17bp_barcode>`；
- `cell_type` 为空、空白或 `NA`：排除；
- 不在 RNA 注释表中的 ALLC：排除；
- ALLC 源目录只读，转换结果写入当前 run 的 `01_cov`。

## 运行与验证

推荐在 `cu03` 提交：

```bash
sbatch --partition=cpu --nodelist=cu03 --cpus-per-task=55 --mem=250G \
  /home/lijia/luozhixiong/IPF_tissue/Scripts/Methscan/run_methscan.sbatch \
  /home/lijia/luozhixiong/IPF_tissue/Results/Methscan/<run_name>
```

拆分提交入口：

```bash
/home/lijia/luozhixiong/IPF_tissue/Scripts/Methscan/submit_methscan_pipeline.sh \
  /home/lijia/luozhixiong/IPF_tissue/Results/Methscan/<run_name>
```

公共作业默认申请 55 CPU/250G；每个 threshold 分支默认申请 18 CPU/80G，由 Slurm 自动调度到可用 CPU 节点。

prepare、filter、smooth 已拆为独立串行作业，各申请 4 CPU/16G；selection+ALLC→COV 使用 cu03 的 55 CPU/250G。这样不会为不支持并行的 QC 工具长期占用转换阶段的大资源。

检查 `stage_status.tsv`、Slurm 日志和每个 `<stage>.resources.json`。最终必须满足：

- `run_summary.json` 存在且 `status=complete`；
- 三个 threshold 分支均有 VMR、四个 matrix、Scanpy h5ad 和图；
- matrix 行 ID 与 filtered header 一致；
- Scanpy embedding 是 filtered cell 的子集；
- 所有资源记录 `return_code=0`。

部分运行或取消的目录不得直接续跑；修复后使用新的输出目录。

## 最近 smoke test 记录

## 正式运行记录

作业 `307549` 于 2026-08-26/27 完成 `CYL_ZCP_full_20260826_final` 全流程。8,949 个 ALLC 中 8,626 个通过 RNA cell-ID 及非空/非 `NA` cell type 筛选；MethSCAn filter 后保留 6,264 个细胞。三个分支分别产生 39,553、80,818、166,618 个 VMR，三个分支均完成 matrix 和 Scanpy 输出，`run_summary.json` 为 `status=complete`。prepare 于 2026-08-26 18:04 完成，最终 summary 于 2026-08-27 02:59 完成。

2026-08-26 在 `/tmp/methscan_smoke_reorg_20260826` 完成了 20 个细胞（CYL/ZCP 各 10）全流程。三个分支均完成：0.01 产生 3,240 个 VMR，0.02 产生 6,784 个 VMR，0.05 产生 22,327 个 VMR；过滤后保留 13 个细胞，最终 summary 为 `complete`。

测试发现并修复了一个确定问题：`05_summarize_run.py` 要求 `02_prepared/cell_id_check.json`，而原 `03_prepare_methscan.py` 未生成该文件。现在 prepare 阶段会比较 manifest、`column_header.txt` 和 `cell_stats.csv` 的 canonical cell ID，并在不一致时失败。
