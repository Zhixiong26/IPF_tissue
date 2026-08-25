# MethSCAn：ALLCools 到 VMR/Scanpy / ALLCools-to-VMR/Scanpy

本阶段从逐细胞 ALLCools ALLC 文件开始，自动完成输入校验、可选 QC 选择、MethSCAn VMR 发现、矩阵构建和 Scanpy 降维聚类。原始输入只读；每次正式运行使用新的 `Results/Methscan/<run_name>/`，日志写入 `Scripts/Methscan/logs/`。

This stage starts from per-cell ALLCools ALLC files and automates input validation, optional QC selection, MethSCAn VMR discovery, matrix construction, and Scanpy embedding/clustering. Inputs remain read-only; every formal run uses a new `Results/Methscan/<run_name>/`, with logs under `Scripts/Methscan/logs/`.

## 已确认输入与环境 / Confirmed inputs and environments

- 当前 ALLC 源 / Current ALLC source: `/home/lijia/jiangyuanpei/methscan/xunyin/IPF_tissue/allcools_5kbin/input_allc`
- 6,554 cells：CYL 3,165；ZCP 3,389；每个 `*.allc.tsv.gz` 均有 `.tbi`。
- 抽查为标准 ALLC 七列，当前只有 `CGN`，所以本流程发现的是 CpG VMR，不提供 mCH/mCCC。
- MethSCAn: `/home/lijia/jiangyuanpei/miniforge3/envs/MethSCAn/bin/methscan`，v1.1.0。
- Scanpy: `/home/lijia/jiangyuanpei/miniforge3/envs/allcools/bin/python`，Scanpy 1.9.3 基线。

The current source has 6,554 indexed ALLCs (3,165 CYL and 3,389 ZCP). Sampled records follow the seven-column ALLC contract and currently contain only `CGN`; this workflow therefore discovers CpG VMRs and cannot recover mCH/mCCC.

## 文件与步骤 / Files and stages

1. `00_methscan_config.sh`: 输入、环境、阈值、线程和可选 TSS BED。 / Inputs, environments, thresholds, threads, and optional TSS BED.
2. `01_prepare_allc_inputs.py`: 接受 ALLC 目录或 `tar(.gz)`；安全解包、格式/索引检查、QC 选择、平衡 smoke 子集、规范化软链接和 manifest。 / Accepts an ALLC directory or tar archive; safely extracts, validates format/index, applies optional QC, builds balanced smoke subsets, canonical links, and a manifest.
3. `02_prepare_methscan.py`: 直接运行 `methscan prepare --input-format allc`，不生成冗余 coverage。 / Calls `methscan prepare --input-format allc` directly without redundant coverage conversion.
4. `preflight_methscan.sbatch`: 在 fat 计算节点抽查 ALLC 路径、索引和两个环境。 / Checks the ALLC path, indices, and both environments from a fat compute node.
5. `run_methscan.sbatch`: 串联 filter → smooth → scan → matrix → Scanpy，并可选执行 TSS profile。 / Chains filter, smooth, scan, matrix, Scanpy, and optional TSS profiling.
6. `03_vmr_scanpy.py`: 缺失率过滤、迭代 PCA 填补、PCA、邻接图、UMAP、Leiden、H5AD/表格/图片。 / Missingness filtering, iterative PCA imputation, PCA, neighbours, UMAP, Leiden, and H5AD/table/figure export.
7. `04_summarize_run.py`: 验证必需输出并生成 `run_summary.json/tsv`。 / Validates required products and writes `run_summary.json/tsv`.

输出批次结构 / Run output layout:

```text
Results/Methscan/<run_name>/
├── 00_manifest/
├── 01_prepared/
├── 02_filtered/
│   └── smoothed/
├── 03_scan/VMRs.bed
├── 04_matrix/
├── 05_scanpy/
├── 06_profile/              # only when TSS BED is configured
├── software_versions.txt
├── run_summary.json
└── run_summary.tsv
```

## QC 接口 / QC interface

QC 结果尚未完成时，不传 QC 表，流程保留全部可用 ALLC，同时仍执行 MethSCAn 的 `min-sites/min-meth/max-meth` 过滤。QC 完成后，把逐细胞主表和明确的布尔列传给提交脚本。当前 ALLC 只有 CGN，若 `pass_final_qc` 因 mCH/mCCC 缺失而为 `NA`，不能盲目使用它；应先确认采用 `pass_CpG`、`pass_mapping` 等列还是新定义的可用指标组合。

When QC is pending, omit the QC table; all available ALLCs enter the workflow and MethSCAn's site/methylation filters still apply. Once QC is complete, pass the per-cell master table and an explicitly reviewed Boolean column. Because current ALLCs contain only CGN, do not use `pass_final_qc` blindly if mCH/mCCC make it `NA`.

## Smoke 与正式运行 / Smoke and formal runs

先提交轻量计算节点预检： / Submit the lightweight compute-node preflight first:

```bash
cd /home/lijia/luozhixiong/IPF_tissue
sbatch Scripts/Methscan/preflight_methscan.sbatch
```

先用平衡的小批次检查实际内存、磁盘和运行时间： / First use a balanced subset to measure memory, storage, and runtime:

```bash
cd /home/lijia/luozhixiong/IPF_tissue
IPF_METHSCAN_MAX_CELLS=20 sbatch \
  Scripts/Methscan/run_methscan.sbatch \
  /home/lijia/luozhixiong/IPF_tissue/Results/Methscan/smoke_YYYYMMDD
```

正式运行（QC 尚未完成）： / Formal run while QC is pending:

```bash
sbatch Scripts/Methscan/run_methscan.sbatch \
  /home/lijia/luozhixiong/IPF_tissue/Results/Methscan/CYL_ZCP_YYYYMMDD
```

QC 完成后的接口示例；提交前必须确认 QC 列： / Example after QC completion; review the QC column before submission:

```bash
sbatch Scripts/Methscan/run_methscan.sbatch \
  /home/lijia/luozhixiong/IPF_tissue/Results/Methscan/CYL_ZCP_QC_YYYYMMDD \
  /path/to/per_cell_qc.tsv \
  pass_final_qc
```

提交脚本拒绝已有结果路径，要求至少 500 GiB 可用空间，并申请 fat 分区 32 CPU、256 GiB、5 天。正式提交前先在计算节点验证 ALLC 源可读。TSS profile 只有在参考基因组和字典序排序的 BED 确认后才通过 `IPF_METHSCAN_TSS_BED` 启用。

The launcher refuses an existing result path, requires at least 500 GiB free, and requests 32 CPUs, 256 GiB, and five days on `fat`. Verify compute-node access to the ALLC source before submission. TSS profiling is enabled through `IPF_METHSCAN_TSS_BED` only after genome build and alphabetically sorted BED validation.
