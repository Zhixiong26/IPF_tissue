# MethSCAn 分析报告 / Analysis report

状态 / Status：初版实现完成，正式运行未提交 / preliminary implementation complete; no formal run submitted

## 已核验事实 / Verified facts

- MethSCAn 1.1.0 的真实 CLI 已核验：`prepare`, `filter`, `smooth`, `scan`, `matrix`, `profile`。
- `prepare` 原生支持 `--input-format allc`，因此主流程不再先转换成 Bismark coverage。
- 正式默认输入已切换为项目内 `Data/ALLCools/` 下的 CYL/ZCP 原始 `allcools.tar.gz` 归档；每次新运行在其 Results 目录中安全解包，原始归档保持只读。归档复制完成后由 intake 重新记录 indexed ALLC 数量、contexts 和索引状态。
- 4-cell intake smoke（CYL/ZCP 各 2）通过；每个文件抽查 100 条，共 400 条 `CGN` 记录，规范 cell ID 和软链接正确。
- canonical RNA 注释 `Results/Scanpy/E_CYL_ZCP_notebook/cell_id_cell_type.tsv` 的 `cell_id` 是 MethSCAn 的唯一 RNA 白名单。与 ALLC cell IDs 的既往交集为 6,317/6,554；任何 MethSCAn 命令前只保留精确匹配且 `cell_type` 不为字面 `NA` 的细胞，未匹配和 `NA` 细胞分别审计并排除。该筛选只执行一次；filter 后仅按甲基化 QC 去除细胞，不再重新匹配或复核 cell type。
- VMR Scanpy 已增加 cell/VMR 缺失率审计表，最终 summary 会严格检查 manifest→prepare→filter→matrix→Scanpy 的细胞数一致性。
- Python 编译、Shell 语法、命令行接口和 skill package 校验通过。Slurm `--test-only` 接受 2-CPU preflight 与 32-CPU full wrapper；测试编号 307509/307510 不是实际任务。尚未运行 MethSCAn prepare 或后续全链路。

The MethSCAn 1.1.0 CLI was verified. Its native ALLC input is used directly. The current source contains 6,554 indexed ALLCs (~19 GB; 3,165 CYL and 3,389 ZCP). A balanced four-cell intake smoke passed, validating 400 CGN records and canonical cell links. The maintained workflow uses `cell_id_cell_type.tsv` as the sole RNA whitelist before prepare; MethSCAn filter subsequently removes cells only by methylation QC and does not repeat cell-type matching. Missingness audit tables and strict cross-stage cell-count consistency checks are implemented. The real full-chain smoke status is recorded below.

## 计算节点与 smoke 审计 / Compute-node and smoke audit

- fat-node preflight `307517`：完成，退出码 `0:0`；20 个真实 ALLC、索引、MethSCAn 1.1.0 和 Scanpy 1.9.3 可用。
- 更新后 preflight `307521`：完成，退出码 `0:0`；open-file limit 131,072、可用空间约 7.0 TB、便携资源记录器可编译，峰值 RSS 约 71 MB。
- smoke `307518` / `smoke_20260825_v2`：在任何数据计算前失败，因为 fat01 没有 `/usr/bin/time`；保留失败根并以 `05_run_with_resources.py` 替代节点依赖。
- smoke `307522` / `smoke_20260825_v3`：在 `prepare` 阶段由用户主动停止，因为启动时误用了参考示例的 60,000 阈值；未据此解释过滤结果。
- smoke `307525` / `smoke_20260825_v4`：使用项目标准 `min-sites=300,000` 运行；最终状态和结果指标在作业完成后填写。

The two fat-node preflights passed. The first full-chain attempt failed safely before data computation because GNU `time` is absent on fat01; a tested Python resource wrapper now replaces that dependency. The next attempt was deliberately cancelled when the covered-CpG threshold discrepancy was identified. Job `307525` is the first valid smoke launched with the project-standard 300,000-site threshold.

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

`min-sites=300,000` 是本项目固定的 covered CpG 技术门槛；用户确认当前唯一甲基化水平门槛为 overall mCG >0.5。MethSCAn 的 `--min-meth` 使用百分数并包含边界，因此实现为 `min-meth=50`；`max-meth=100` 不施加有效上限。20-cell smoke 的 overall mCG 为 0.6946–0.7726，按新门槛回看仍为 20/20 通过，但该 smoke 的作业记录仍保留其实际运行参数 20–85%。 / `min-sites=300,000` remains the technical covered-CpG threshold. The user confirmed that the only methylation-level gate is overall mCG >0.5. MethSCAn uses percentages with an inclusive minimum, implemented as `min-meth=50`; `max-meth=100` is nonrestrictive. The 20-cell smoke spans 0.6946–0.7726 overall mCG and retrospectively retains all 20 under the new threshold, while its audit record continues to report the actual 20–85% parameters used at runtime.

## 未决事项 / Open items

- 不使用已删除的 QC 主表。Scanpy `cell_type=NA` 是明确排除条件，不得作为可进入 MethSCAn 的细胞。
- 参考基因组尚未确认，TSS BED/profile 默认关闭。
- 正式分析前需要计算节点读取 ALLC 源的 preflight、空间评估和 20-cell smoke。
- 还没有 VMR 数量、过滤保留率、缺失率、PCA/UMAP/Leiden 或生物学注释结果。

The genome build/TSS BED remain unresolved. The user explicitly chose not to rerun the 20-cell smoke after adding the RNA-first gate; any later user-authorized run applies that gate directly. No VMR, retention, missingness, embedding, clustering, or biological annotation result exists yet.
