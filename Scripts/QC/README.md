# Per-cell QC / 逐细胞质控

本阶段按 FASTQ、BAM+CB、ALLC 三层恢复指标，最后生成逐细胞主表和分样本汇总。所有 PASS/FAIL 标志均保留；脚本不会删除细胞。

This stage recovers metrics independently from FASTQ, BAM+CB, and ALLC, then produces a per-cell master table and per-sample summary. Every PASS/FAIL flag is retained; no script deletes cells.

## 环境 / Environment

脚本只使用 Python 标准库。BAM 扫描调用已验证的 Samtools 1.22：

The scripts use only the Python standard library. BAM scanning calls the verified Samtools 1.22 executable:

```bash
PYTHON=/home/lijia/jiangyuanpei/miniforge3/envs/allcools/bin/python
SAMTOOLS=/home/lijia/jiangyuanpei/miniforge3/envs/allcools/bin/samtools
```

结果写入独立的 `Results/QC/<run_name>/`，日志写入 `Scripts/QC/logs/`。不要复用已有 run root。

Write results to a distinct `Results/QC/<run_name>/` and logs to `Scripts/QC/logs/`. Do not reuse an existing run root.

## 输入约束 / Input contract

- FASTQ 必须成对，R1/R2 read name 一致。本项目已确认 `R1[0:17]` 是 cell barcode，`R1[17:29]` 是 UMI；header 末尾 8 bp 是 sample index。 / FASTQs must be paired with matching R1/R2 read names. This project has confirmed that `R1[0:17]` is the cell barcode and `R1[17:29]` is the UMI; the 8-bp header suffix is the sample index.
- BAM 必须含 `CB:Z:<barcode>`；统计时排除 unmapped、secondary、supplementary records，`mapped_pairs` 只数 read1 (`FLAG & 64`)。 / BAMs must contain `CB:Z:<barcode>`; unmapped, secondary, and supplementary records are excluded, and `mapped_pairs` counts read1 only (`FLAG & 64`).
- ALLC manifest 是带表头的 TSV，列为 `sample_id`, `cell_id`, `allc_path`，可选 `barcode`。 / The ALLC manifest is a headered TSV with `sample_id`, `cell_id`, and `allc_path`, optionally `barcode`.
- 现有项目 ALLC 由六列 Bismark CpG coverage 转换，只含 `CGN`。它们可提供 mCG/CpG coverage，但不能提供 mCH/mCCC；脚本将后二者写为 `NA`。 / Current project ALLCs were converted from six-column Bismark CpG coverage and contain only `CGN`. They provide mCG/CpG coverage but not mCH/mCCC; the latter are written as `NA`.

## 执行顺序 / Execution order

```bash
$PYTHON Scripts/QC/01_FASTQ_per_cell_reads.py --help
$PYTHON Scripts/QC/02_BAM_per_cell_mapping.py --help
$PYTHON Scripts/QC/03_merge_FASTQ_BAM.py --help
$PYTHON Scripts/QC/04_ALLC_methylation_QC.py --help
$PYTHON Scripts/QC/05_merge_all_QC.py --help
$PYTHON Scripts/QC/06_QC_summary.py --help
```

01 与 02 是数据层编号。为避免把测序错误产生的 barcode 当作 cell，实际执行时先运行 02 获得 BAM barcode 表，再把它作为 01 的 whitelist。02 可用 `--bam-dir` 自动发现一个样本的全部 Bismark BAM。

Steps 01 and 02 are layer labels. In execution, run 02 first to obtain the BAM barcode table, then use it as the step-01 whitelist so sequencing-error barcodes are not treated as cells. Step 02 can discover all Bismark BAMs for one sample through `--bam-dir`.

```bash
$PYTHON Scripts/QC/02_BAM_per_cell_mapping.py \
  --sample-id CYL \
  --bam-dir Data/Bam/25100718_CYL_Met/bismark \
  --samtools "$SAMTOOLS" \
  --output Results/QC/RUN/02_CYL_bam.tsv

$PYTHON Scripts/QC/01_FASTQ_per_cell_reads.py \
  --sample-id CYL \
  --r1 Data/Raw_fastq/Met/25100718_CYL_Met/1/25100718_CYL_Met_S01_L000_R1_001.fastq.gz \
  --r2 Data/Raw_fastq/Met/25100718_CYL_Met/1/25100718_CYL_Met_S01_L000_R2_001.fastq.gz \
  --r1 Data/Raw_fastq/Met/25100718_CYL_Met/2/25100718_CYL_Met_S01_L003_R1_001.fastq.gz \
  --r2 Data/Raw_fastq/Met/25100718_CYL_Met/2/25100718_CYL_Met_S01_L003_R2_001.fastq.gz \
  --whitelist Results/QC/RUN/02_CYL_bam.tsv \
  --output Results/QC/RUN/01_CYL_fastq.tsv
```

01 的项目默认参数已经是 `--barcode-source r1 --barcode-start 0 --barcode-length 17`。白名单既可是一行一个 barcode 的文本文件，也可直接使用带 `sample_id/barcode` 列的 02 输出；脚本会自动只选当前 sample。

The project defaults for step 01 are now `--barcode-source r1 --barcode-start 0 --barcode-length 17`. The whitelist may be a one-barcode-per-line file or the step-02 output containing `sample_id/barcode`; the script automatically selects the current sample.

## 并行策略 / Parallel strategy

- `01_FASTQ_per_cell_reads.py --workers 2`：并行扫描一个样本的两个 FASTQ partitions。 / Scans the two FASTQ partitions of one sample concurrently.
- `02_BAM_per_cell_mapping.py --workers 16`：每个纳入分析的样本有 156 个 prefix BAM（CYL 约 0.54 TiB，ZCP 约 0.63 TiB），默认使用 16 个并发 Samtools streams；两个 array tasks 同时运行时共 32 workers。 / Each included sample has 156 prefix BAMs (about 0.54 TiB for CYL and 0.63 TiB for ZCP). The default is 16 concurrent Samtools streams, or 32 workers when both array tasks run together.
- `04_ALLC_methylation_QC.py --workers 16`：按逐细胞 ALLC 文件并行；需要时可在独立 benchmark 后提升到 24/32。 / Processes per-cell ALLC files concurrently; it may be raised to 24/32 after an independent benchmark when needed.
- `03/05/06` 只合并小型 TSV，保持串行以避免额外调度开销。 / Steps 03/05/06 only merge small TSVs and remain serial to avoid scheduler overhead.

三个并行脚本均支持 `--progress-interval SECONDS`，默认及正式 `07`/`08` Slurm 入口均设为每 300 秒（5 分钟）报告已完成任务数和运行时间。heartbeat 只汇总已完成的 BAM、FASTQ partition 或 ALLC 文件，不额外扫描输入数据。

All three parallel scripts support `--progress-interval SECONDS`. Both the default and the formal step-07/08 Slurm entry points report completed tasks and elapsed time every 300 seconds (5 minutes). The heartbeat only summarizes completed BAMs, FASTQ partitions, or ALLC files and does not rescan input data.

为避免 BAM 和 FASTQ 使用相同资源造成浪费，集群入口拆为两个相互依赖的两样本 arrays，仅运行当前分析纳入的 CYL 和 ZCP。两个 array 均提交到 `fat` partition：`07` 每个样本申请 16 CPUs/32 GiB，`08` 每个样本只申请 2 CPUs/8 GiB。两个脚本都要求显式指定具有足够空间的新 run root。

To avoid wasting the same allocation on dissimilar BAM and FASTQ stages, the cluster entry point is split into two dependent two-sample arrays covering only CYL and ZCP, the samples included in the current analysis. Both arrays use the `fat` partition. Step 07 requests 16 CPUs/32 GiB per sample and step 08 requests only 2 CPUs/8 GiB per sample. Both require an explicit new run root with sufficient space.

```bash
OUTPUT_ROOT=/home/lijia/luozhixiong/IPF_tissue/Data/QC/qc_counts_RUN_NAME
BAM_JOB=$(sbatch --parsable Scripts/QC/07_run_parallel_bam.sbatch "$OUTPUT_ROOT")
sbatch --dependency="aftercorr:${BAM_JOB}" Scripts/QC/08_run_parallel_fastq.sbatch "$OUTPUT_ROOT"
```

脚本拒绝覆盖已有的 `01_fastq.tsv` 或 `02_bam.tsv`，并要求 `OUTPUT_ROOT` 位于项目 `Data/` 下。当前正式输出根为 `/home/lijia/luozhixiong/IPF_tissue/Data/QC/full_counts_20260824`。

The scripts refuse to overwrite existing `01_fastq.tsv` or `02_bam.tsv` and require `OUTPUT_ROOT` to be under the project `Data/` directory. The current formal output root is `/home/lijia/luozhixiong/IPF_tissue/Data/QC/full_counts_20260824`.

这里没有把每个样本设为 56 workers：两样本 BAM 同时使用 32 个并发读取流，瓶颈预期是共享 `/mnt/data04` 的吞吐而不是 CPU。若监测显示 I/O wait 较低且 CPU 仍空闲，再把 `--cpus-per-task` 和 `--workers` 一起提升到 24；不要只增加 workers 而不增加 Slurm CPU allocation。

The design intentionally does not use 56 workers per sample: the two BAM samples already create 32 concurrent read streams, and shared `/mnt/data04` throughput is expected to bottleneck before CPU. If monitoring shows low I/O wait and spare CPU, raise both `--cpus-per-task` and `--workers` to 24; do not increase workers without increasing the Slurm CPU allocation.

## 阈值 / Thresholds

`05` 默认使用严格阈值：`500,000 < final_reads < 10,000,000`、`mapping_rate > 0.50`、`mCG > 0.50`、`mCH < 0.20`、`mCCC < 0.05`。coverage secondary QC 使用闭区间 `300,000 <= n_CpG_covered <= 1,200,000`。所有阈值均可通过命令行覆盖。

`05` defaults to strict thresholds: `500,000 < final_reads < 10,000,000`, `mapping_rate > 0.50`, `mCG > 0.50`, `mCH < 0.20`, and `mCCC < 0.05`. The secondary coverage QC uses the inclusive interval `300,000 <= n_CpG_covered <= 1,200,000`. Every threshold is configurable on the command line.

`final_reads` 当前定义为 BAM 中带 CB 的 primary mapped R1+R2 records，并同时保留 `mapped_pairs` 与 `mapped_reads`。论文定义确认后可调整，不必重扫 BAM。

`final_reads` is currently defined as primary mapped R1+R2 BAM records carrying CB, while `mapped_pairs` and `mapped_reads` are retained. The definition can be revised after checking the paper without rescanning BAM.

`03` 的 `--fastq` 与 `--bam` 可以各自重复两次，从而直接合并 CYL、ZCP 的逐样本结果。`n_CpG_covered` 的精确定义是 coverage > 0 的 CG-context ALLC 行数；若将来提供 strand-split ALLC，需要先确认论文所称 CpG site 是否要求合并互补链。

The `--fastq` and `--bam` arguments to `03` can each be repeated twice to combine CYL and ZCP outputs directly. `n_CpG_covered` is defined exactly as the number of CG-context ALLC rows with coverage > 0; for future strand-split ALLCs, confirm whether the paper's CpG-site definition requires complementary-strand merging.
