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
07_methdiff_celltype.py
08_hypo_dmr_heatmaps.py
```

主入口 `run_methscan.sbatch` 已同步更新引用，不再调用不存在的 `01_prepare_allc_inputs.py`，也不再使用旧的 `link_path`/`00_manifest` intake 阶段。

新增 `run_methscan_common.sbatch`、`run_methscan_qc_stage.sbatch`、`run_methscan_branch.sbatch`、`run_methscan_summary.sbatch` 和 `submit_methscan_pipeline.sh`。后续正式运行只需执行一次 `submit_methscan_pipeline.sh`：它自动提交公共阶段、三个独立 QC 阶段、三个 threshold 分支和 summary，并设置 `afterok` 前后依赖。prepare/filter/smooth 各申请 4 CPU/16G。

## 当前流程

```text
ALLC → Scanpy whitelist → ALLC→cov → prepare → filter → smooth
→ 样本内 pairwise cell-type meth-diff
→ 可选 raw-p fallback → Top200 hypo-DMR → 单细胞×DMR矩阵
→ raw mean ratio / z-score / 压缩 colorbar z-score 热图
→ (scan → matrix → Scanpy) × {0.01, 0.02, 0.05} → summary

独立可选分支：
smooth → 合并 CYL+ZCP 同名 cell type → 91 个 pooled pairwise meth-diff
```

pairwise meth-diff 使用 smooth 后 `03_filtered` 的 chr1–22/X/Y hard-link view（不复制数据）。主流程 group 文件按 sample 限定细胞，在 `07_methdiff/samples/<sample>/comparisons/` 下写入结果；独立 pooled 分支把 CYL/ZCP 同名 cell type 合并后写入 `pooled/07_methdiff/samples/CYL_ZCP/`。分样本热图下游严格采用参考流程的 raw p `<0.01`、mean ratio 绝对差 `≥0.25`、每细胞类型 Top200 hypo-DMR 口径；矩阵来自 COV 内 unique CpG ratio 的等权平均。三套热图分别为原始 mean ratio、DMR-wise z-score，以及仅把 colorbar 压缩到 `[-1,1]` 的未改值 z-score，并显示行列 cell-type 标签与分组边界。

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

prepare、filter、smooth 已拆为独立串行作业，各申请 4 CPU/16G。pairwise meth-diff 默认以 `2 comparisons × 4 threads` 滚动执行，申请 8 CPU/16G；`07_methdiff_celltype.py` 强制检查 `jobs × threads ≤ SLURM_CPUS_PER_TASK`，各比较使用独立进程、独立输出且 MethSCAn 固定随机种子，因而并发不改变单比较统计口径。hypo-DMR 矩阵/绘图默认使用 8 workers、8 CPU/48G。两者允许在 `cpu,fat` 分区间择空调度，并可通过 `IPF_METHSCAN_METHDIFF_JOBS`、`IPF_METHSCAN_METHDIFF_THREADS`、`IPF_METHSCAN_HYPO_MATRIX_WORKERS` 或 `sbatch` 资源参数覆盖。selection+ALLC→COV 使用 cu03 的 55 CPU/250G。

检查 `stage_status.tsv`、Slurm 日志和每个 `<stage>.resources.json`。最终必须满足：

- `run_summary.json` 存在且 `status=complete`；
- `07_methdiff/pairwise_summary.json` 不含硬失败；已知 FDR 除零比较由 `08_rawp_fallback` 定向补跑；
- `hypo_dmr_heatmap_summary.json` 存在且 `status=complete`；
- 三个 threshold 分支均有 VMR、四个 matrix、Scanpy h5ad 和图；
- matrix 行 ID 与 filtered header 一致；
- Scanpy embedding 是 filtered cell 的子集；
- 所有资源记录 `return_code=0`。

部分运行或取消的目录不得直接续跑；修复后使用新的输出目录。

## 正式运行记录

作业 `307549` 于 2026-08-26/27 完成 `CYL_ZCP_full_20260826_final` 全流程。8,949 个 ALLC 中 8,626 个通过 RNA cell-ID 及非空/非 `NA` cell type 筛选；MethSCAn filter 后保留 6,264 个细胞。三个分支分别产生 39,553、80,818、166,618 个 VMR，三个分支均完成 matrix 和 Scanpy 输出，`run_summary.json` 为 `status=complete`。prepare 于 2026-08-26 18:04 完成，最终 summary 于 2026-08-27 02:59 完成。

分样本 hypo-DMR 下游已完成：CYL 为 2,409 个合并 Top200 hypo-DMR × 2,919 个细胞，ZCP 为 2,502 个合并 Top200 hypo-DMR × 3,345 个细胞；两者均生成带 cell-type 行列标签的 mean-ratio、z-score 和压缩 colorbar z-score 热图。

pooled 作业 `309188` 于 2026-09-16 完成，耗时 5:23:06。它将 CYL/ZCP 的同名 RNA cell type 合并，在全部 6,264 个过滤后细胞中完成 14 种细胞类型的 91/91 个两两比较；0 失败、0 fallback、0 空结果，共 6,437,767 行 DMR。pooled 阶段只负责 MethSCAn DMR，不生成 Top200、单细胞矩阵或热图；其全量 unique-DMR 后续由 `Scripts/Methylvi/vmr_dmr/` 消费。
