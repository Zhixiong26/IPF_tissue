# FastQ 分析报告 / Analysis report

状态 / Status：2026-08-24 已通过轻量级 Raw FASTQ 验证 / lightweight validation passed

## 输入 / Inputs

- `Data/Raw_fastq/E/`：3 个链接批次 / 3 linked batches.
- `Data/Raw_fastq/Met/`：3 个链接批次 / 3 linked batches.
- 每个批次包含子集 `1/`、`2/` 及配对的 R1/R2 文件。 / Each batch contains subsets `1/` and `2/` with paired R1/R2 files.

## 已验证结果 / Verified results

- 6 个有效源链接。 / 6 valid source links.
- 24 个可读非空 FASTQ，组成 12 对完整 R1/R2。 / 24 readable, non-empty FASTQs forming 12 complete R1/R2 pairs.
- 所有文件首条记录的 gzip 解码和 FASTQ 格式均通过。 / First-record gzip decoding and FASTQ syntax passed for every file.
- 每对文件的首条读名一致。 / First-record read names matched within every pair.
- 总压缩大小 / Total compressed size：1,663,551,234,393 bytes.

机器可读结果位于 `Results/FastQ/raw_fastq_validation/`。 / Machine-readable outputs are stored there.

## 限制与待确认事项 / Limitations and pending decisions

- 未验证全流 gzip 完整性和全文件校验和，因为这需要顺序读取约 1.66 TB。 / Full-stream gzip integrity and whole-file checksums were not tested because that would sequentially read about 1.66 TB.
- Bismark 设计前仍需确认 `E`、`Met`、子集 `1/2` 的含义、建库方案和参考基因组。 / The meanings of `E`, `Met`, and subsets `1/2`, the library protocol, and reference genome remain to be confirmed.

## 已暂缓合并 / Deferred merge

- 计划输出 / Planned output: `Results/FastQ/E/25100718_CYL_E/25100718_CYL_E_R1_001.fastq.gz`.
- 两个 R1 输入合计 21,723,176,838 bytes。 / The two R1 inputs total 21,723,176,838 bytes.
- 2026-08-24 检查时目标文件系统仅余约 6.2 GB，因此未开始合并，也未产生残缺文件。 / At the 2026-08-24 check, the target filesystem had only about 6.2 GB free, so the merge was not started and no partial output was created.
- 使用者已明确选择当前项目从 `Data/Bam` 开始，FASTQ 合并和处理暂时跳过。 / The user explicitly chose to start this project from `Data/Bam`; FASTQ merging and processing are temporarily deferred.
