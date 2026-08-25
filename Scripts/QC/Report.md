# QC 分析报告 / QC analysis report

状态 / Status: barcode 规则及脚本已验证，待 CYL、ZCP 两样本完整计数；full-context ALLC 独立待办 / barcode rule and scripts verified; full counting for CYL and ZCP pending, with full-context ALLC tracked separately

## 已完成 / Completed

- 建立 `01`–`06` 六步逐细胞 QC 脚本，分别统计 FASTQ、BAM、ALLC、合并、阈值标记及分样本汇总。 / Implemented six per-cell QC steps for FASTQ, BAM, ALLC, merging, threshold flags, and per-sample summaries.
- BAM 代表记录确认含 17-bp `CB:Z`，paired-end Bismark records 可通过 read1 flag 统计 mapped pairs。 / Representative BAM records contain a 17-bp `CB:Z`; paired-end Bismark mapped pairs can be counted through the read1 flag.
- ALLC 代表文件为七列格式，但 context 仅观察到 `CGN`。 / Representative ALLC files have seven columns, but only `CGN` context was observed.
- QC 标记缺失值使用 `NA`，不把缺少 mCH/mCCC 的 CpG-only ALLC 当作 0 methylation。 / Missing QC values are represented as `NA`; CpG-only ALLCs lacking mCH/mCCC are not treated as zero methylation.
- 已确认 Met R1 前 17 bp 是 CB、随后 12 bp 是 UMI；FASTQ header 的 8-bp 尾序列是 sample index。 / Confirmed that the first 17 bases of Met R1 are CB and the next 12 are UMI; the 8-bp FASTQ-header suffix is the sample index.

## 当前限制 / Current limitations

- 现有 CpG-only ALLC 无法计算 mCH/mCCC。完成 primary QC 需从保留非-CG context 的 methylation calls 重新生成 ALLC，或提供等价逐细胞统计。 / Existing CpG-only ALLCs cannot provide mCH/mCCC. Completing primary QC requires ALLCs rebuilt from methylation calls retaining non-CG contexts, or equivalent per-cell metrics.
- 尚未批量扫描当前纳入分析的 CYL、ZCP 大型 BAM/FASTQ；当前阶段只实现并进行小规模验证。 / The large BAM/FASTQ inputs for CYL and ZCP, the samples currently included in the analysis, have not yet been scanned in full; this stage implements the workflow and performs small-scale validation only.

## 验证 / Verification

- `Results/QC/smoke_20260824/` 的两细胞端到端测试得到预期结果：一个通过 primary/final QC，另一个因 mCG、mCH、mCCC 失败。 / The two-cell end-to-end test under `Results/QC/smoke_20260824/` produced the expected result: one primary/final pass and one mCG/mCH/mCCC failure.
- CG-only fixture 正确将 mCH/mCCC 写为 `NA`，并标记 `CG_only_CH_CCC_unavailable`。 / The CG-only fixture correctly writes mCH/mCCC as `NA` and marks `CG_only_CH_CCC_unavailable`.
- 真实 CYL BAM 前 1,000 条 primary mapped records：CB 缺失 0，read1 500，read2 500，观察到 108 个 barcode。 / In the first 1,000 primary mapped records of a real CYL BAM: 0 lacked CB, 500 were read1, 500 were read2, and 108 barcodes were observed.
- 三个样本分别交叉匹配 50 个 raw/BAM read names：150/150 的 raw R1 前 17 bp 等于 BAM CB，150/150 的随后 12 bp 等于 BAM UR/UMI。 / Fifty raw/BAM read names were cross-matched per sample: 150/150 raw R1 17-bp prefixes equalled BAM CB and 150/150 following 12-bp segments equalled BAM UR/UMI.
- CYL partition 2 前 100,000 对 reads 使用 partial BAM whitelist 的测试成功：barcode 结构异常 0，对 108 个 whitelist cells 中的 80 个分配到 1,196 pairs。 / A 100,000-pair CYL partition-2 test with a partial BAM whitelist succeeded: 0 structurally invalid barcodes and 1,196 pairs assigned across 80 of 108 whitelist cells.
- 扫描密集型步骤已改为有界多进程：01 按 FASTQ partition、02 按 prefix BAM、04 按 cell ALLC 并行；03/05/06 保持串行。 / Scan-heavy steps now use bounded multiprocessing: step 01 by FASTQ partition, step 02 by prefix BAM, and step 04 by cell ALLC; steps 03/05/06 remain serial.
- 01/02/04 增加可配置 heartbeat；正式 07/08 作业每 300 秒（5 分钟）向 stdout 报告完成数/总数和运行时间，便于区分长时间扫描与任务停滞。 / Steps 01/02/04 now provide a configurable heartbeat; formal step-07/08 jobs report completed/total tasks and elapsed time to stdout every 300 seconds (5 minutes), helping distinguish long scans from stalled jobs.
- 使用 2 个 FASTQ partitions、2 个各含 1,000 records 的 BAM fixtures、2 个 ALLC cells 比较 `--workers 1` 与 `--workers 2`，三层 TSV 均逐字节一致。 / Comparing `--workers 1` versus `--workers 2` on two FASTQ partitions, two 1,000-record BAM fixtures, and two ALLC cells produced byte-identical TSVs at all three layers.
- 资源复核显示每个样本有 156 个 prefix BAM（CYL 0.537 TiB、LC 0.498 TiB、ZCP 0.632 TiB），但每个样本只有 2 个 R1 FASTQ partitions。BAM 默认提升为 16 workers，ALLC 提升为 16 workers；FASTQ 保持 2。 / Resource review found 156 prefix BAMs per sample (CYL 0.537 TiB, LC 0.498 TiB, ZCP 0.632 TiB) but only two R1 FASTQ partitions per sample. BAM and ALLC defaults were raised to 16 workers; FASTQ remains at two.
- 原组合脚本已替换为未提交的 `07_run_parallel_bam.sbatch`（`fat` partition，2 × 16 CPUs/32 GiB array）和依赖的 `08_run_parallel_fastq.sbatch`（`fat` partition，2 × 2 CPUs/8 GiB array）；正式范围仅含 CYL 和 ZCP，避免 FASTQ 阶段空占 14 CPUs/样本。 / The former combined script was replaced by the non-submitted `07_run_parallel_bam.sbatch` (`fat` partition, 2 × 16 CPUs/32 GiB array) and dependent `08_run_parallel_fastq.sbatch` (`fat` partition, 2 × 2 CPUs/8 GiB array). The formal scope contains only CYL and ZCP and avoids 14 idle CPUs per sample during FASTQ scanning.
- 首次正式提交的 BAM array `307475` 在 `fat01` 启动后立即失败，因为输出根位于该节点不可访问的 `/mnt/data04`；依赖的 FASTQ array `307476` 因 `DependencyNeverSatisfied` 已取消。未产生 BAM/FASTQ 计数结果。 / The first formal BAM array submission, `307475`, failed immediately on `fat01` because its output root was under `/mnt/data04`, which is unavailable on that node. The dependent FASTQ array `307476` was cancelled after entering `DependencyNeverSatisfied`. No BAM/FASTQ count result was produced.
- 输出根已改为项目内的 `/home/lijia/luozhixiong/IPF_tissue/Data/QC/full_counts_20260824`，07/08 脚本现在拒绝项目 `Data/` 之外的输出路径。 / The output root is now `/home/lijia/luozhixiong/IPF_tissue/Data/QC/full_counts_20260824`, and step-07/08 now reject output paths outside the project `Data/` directory.
- 当前输入仍阻塞重新提交：`Data/Bam` 和 `Data/Raw_fastq/Met/*` 是指向 `/mnt/data04` 的绝对软连接；`fat01` 只读探测确认这些连接的目标不存在，BAM/FASTQ 均不可读。必须先把真实输入放到 fat01 可访问的项目 `Data/` 路径，或将软连接改到 fat01 可访问的目标。 / Input access still blocks resubmission: `Data/Bam` and `Data/Raw_fastq/Met/*` are absolute symlinks into `/mnt/data04`; a read-only probe on `fat01` confirmed that their targets do not exist there and neither BAM nor FASTQ is readable. The real inputs must first be placed under a project `Data/` path visible to `fat01`, or the symlinks must be changed to targets visible to `fat01`.
- 2026-08-25 的一次提交 `307500` 错误沿用了旧的 `Results/QC/full_counts_20260825` 输出根；当前 launcher 的 Data-prefix guard 使 CYL/ZCP 两个 BAM tasks 在 0 秒内以 exit 3 安全退出，未读取数据。依赖 FASTQ array `307501` 随后取消。 / A 2026-08-25 submission, `307500`, incorrectly reused the old `Results/QC/full_counts_20260825` output root. The current launcher's Data-prefix guard caused both CYL/ZCP BAM tasks to exit safely with code 3 in zero seconds without reading data. Dependent FASTQ array `307501` was then cancelled.
- 07/08 已改为直接读取通过复制校验的 `Data/.fat_stage_20260824/Bam/` 与 `Raw_fastq/Met/`。fat01 上 1 CPU/1 GiB 只读预检确认 CYL/ZCP 各 156 BAM、2 R1、2 R2、环境可执行且 `Data/QC` 可写。 / Steps 07/08 now read directly from the copy-verified `Data/.fat_stage_20260824/Bam/` and `Raw_fastq/Met/`. A 1-CPU/1-GiB read-only preflight on fat01 confirmed 156 BAMs, 2 R1s, and 2 R2s for each CYL/ZCP sample, executable environments, and writable `Data/QC`.
- 正式重提为 BAM array `307505` 和依赖 FASTQ array `307506`，输出根 `Data/QC/full_counts_20260825`。首次核对时 `307505_0/1` 均在 fat01 RUNNING，`307506_0/1` 因正常 dependency PENDING；输出目录和 5 分钟 heartbeat 已建立。此状态仅表示运行开始，不表示 QC 完成。 / Formal resubmission created BAM array `307505` and dependent FASTQ array `307506`, targeting `Data/QC/full_counts_20260825`. At first verification, `307505_0/1` were RUNNING on fat01 and `307506_0/1` were normally PENDING on dependency; output directories and five-minute heartbeat logs were established. This indicates a valid start, not completed QC.

## 不变量 / Invariants

脚本仅生成表格和指标，不删除任何细胞；各项失败可以重叠，因此 summary 中 fail 列不能相加解释为总删除数。

The scripts only generate tables and metrics and never delete cells. Failure categories may overlap, so summary failure columns must not be summed as a deletion total.
