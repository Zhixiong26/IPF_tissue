# MethSCAn 分析报告 / Analysis report

状态 / Status：初版实现完成，正式运行未提交 / preliminary implementation complete; no formal run submitted

## 已核验事实 / Verified facts

- MethSCAn 1.1.0 的真实 CLI 已核验：`prepare`, `filter`, `smooth`, `scan`, `matrix`, `profile`。
- `prepare` 原生支持 `--input-format allc`，因此主流程不再先转换成 Bismark coverage。
- 正式默认输入已切换为项目内 `Data/ALLCools/` 下的 CYL/ZCP 原始 `allcools.tar.gz` 归档；每次新运行在其 Results 目录中安全解包，原始归档保持只读。归档复制完成后由 intake 重新记录 indexed ALLC 数量、contexts 和索引状态。
- 旧输入曾用于开发阶段验证；其结果不属于当前正式分析流程。
- canonical RNA 注释 `Results/Scanpy/E_CYL_ZCP_notebook/cell_id_cell_type.tsv` 的 `cell_id` 是 MethSCAn 的唯一 RNA 白名单。与 ALLC cell IDs 的既往交集为 6,317/6,554；任何 MethSCAn 命令前只保留精确匹配且 `cell_type` 不为字面 `NA` 的细胞，未匹配和 `NA` 细胞分别审计并排除。该筛选只执行一次；filter 后仅按甲基化 QC 去除细胞，不再重新匹配或复核 cell type。
- VMR Scanpy 已增加 cell/VMR 缺失率审计表，最终 summary 会严格检查 manifest→prepare→filter→matrix→Scanpy 的细胞数一致性。
- Python 编译、Shell 语法、命令行接口和 skill package 校验通过。正式运行使用 32-CPU full wrapper。

The MethSCAn 1.1.0 CLI was verified and its native ALLC input is used directly. The current formal source is the extracted CYL/ZCP data under `Data/ALLCools`. The maintained workflow uses `cell_id_cell_type.tsv` as the sole RNA whitelist before prepare; MethSCAn filter subsequently removes cells only by methylation QC and does not repeat cell-type matching. Missingness audit tables and strict cross-stage cell-count consistency checks are implemented.

## 计算节点与正式运行 / Compute-node and formal run

正式运行使用 fat 分区、32 CPU、256 GiB 和最长 5 天；资源由 `05_run_with_resources.py` 记录。

## 初版固定参数 / Preliminary fixed parameters

| Parameter | Value |
|---|---:|
| MethSCAn min sites / covered CpG sites | 300,000 |
| Overall mCG minimum | 50% (MethSCAn inclusive minimum) |
| Overall mCG maximum | 100% (nonrestrictive domain ceiling) |
| Smooth bandwidth | 1,000 bp |
| Scan bandwidth / step | 2,000 / 100 bp |
| Variable-window fraction | 0.02 |
| Minimum cells per VMR | 6 |
| Scan/matrix threads | 32 |
| VMR observed-cell fraction | 0.05 |
| Minimum observed VMRs/cell | 100 |
| PCA / neighbours / Leiden | 30 / 15 / 0.8 |
| Random seed | 20260825 |

`min-sites=300,000` 是本项目固定的 covered CpG 技术门槛；用户确认当前唯一甲基化水平门槛为 overall mCG >0.5。MethSCAn 的 `--min-meth` 使用百分数并包含边界，因此实现为 `min-meth=50`；`max-meth=100` 不施加有效上限。 / `min-sites=300,000` remains the technical covered-CpG threshold. The user confirmed that the only methylation-level gate is overall mCG >0.5. MethSCAn uses percentages with an inclusive minimum, implemented as `min-meth=50`; `max-meth=100` is nonrestrictive.

## 未决事项 / Open items

- 不使用已删除的 QC 主表。Scanpy `cell_type=NA` 是明确排除条件，不得作为可进入 MethSCAn 的细胞。
- TSS profile 使用用户提供的 `Supplementary/human_hg38_TSS.bed`。原文件共 42,024 行、含第 6 列 strand，但采用自然染色体顺序；为满足 MethSCAn 的字典序要求，默认指向仅排序的 `Supplementary/human_hg38_TSS.methscan.bed`（SHA-256 `90476edd228b8c499af62e52f2fdbc18954e56f90c9d197eb1dccbfca17de5c5`）。原文件未修改。
- 原始 CYL/ZCP 归档在样本目录中解压；正式分析直接使用解压后的 ALLC。
- 当前归档入口、Scanpy 白名单和当前脚本将在正式运行中完成端到端验证；历史开发结果不可复用于当前流程。

TSS profiling uses the user-supplied `Supplementary/human_hg38_TSS.bed`. The original has 42,024 records and strand in column 6 but uses natural chromosome order; the workflow defaults to the sort-only `Supplementary/human_hg38_TSS.methscan.bed` required by MethSCAn (SHA-256 `90476edd228b8c499af62e52f2fdbc18954e56f90c9d197eb1dccbfca17de5c5`). The source file remains unchanged. The extracted CYL/ZCP inputs are used directly from their sample directories. Historical development results must not be reused as current evidence.
