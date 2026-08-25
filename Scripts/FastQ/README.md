# FastQ

本阶段在任何剪切或比对前，对项目的 Raw FASTQ 输入进行清点和验证。

This stage inventories and validates the project's Raw FASTQ inputs before any trimming or alignment.

## 输入与输出 / Inputs and outputs

- 输入 / Input：由使用者提供或明确确认的 `Data/Raw_fastq/`。
- 结果 / Results：`Results/FastQ/`。
- 日志 / Logs：引入调度或长时检查后写入 `Scripts/FastQ/logs/`。

## 执行入口 / Entry points

- `01_validate_raw_fastq.py`：验证源链接、文件可读性和大小、R1/R2 完整性、首条记录的 gzip/FASTQ 格式及首条读名一致性；不执行全流解压。 / Validates source links, readability, size, R1/R2 completeness, first-record gzip/FASTQ syntax, and first read-name agreement; it does not perform full-stream decompression.
- skill 辅助脚本 `.codex/skills/ipf-methylome-analysis/scripts/inspect_fastq_layout.sh` 生成轻量级清单。 / The skill helper produces a lightweight inventory.

使用兼容服务器 Python 3.6 的入口运行验证。 / Run the validator with the server Python 3.6-compatible entry point:

```bash
python3 Scripts/FastQ/01_validate_raw_fastq.py
```
