# MethSCAn 分析报告 / Analysis report

状态 / Status：初版实现完成，正式运行未提交 / preliminary implementation complete; no formal run submitted

## 已核验事实 / Verified facts

- MethSCAn 1.1.0 的真实 CLI 已核验：`prepare`, `filter`, `smooth`, `scan`, `matrix`, `profile`。
- `prepare` 原生支持 `--input-format allc`，因此主流程不再先转换成 Bismark coverage。
- 当前 ALLCools 输入共 6,554 个 indexed ALLC，约 19 GB；CYL 3,165，ZCP 3,389。
- 4-cell intake smoke（CYL/ZCP 各 2）通过；每个文件抽查 100 条，共 400 条 `CGN` 记录，规范 cell ID 和软链接正确。
- Python 编译、Shell 语法、命令行接口和 skill package 校验通过。Slurm `--test-only` 接受 2-CPU preflight 与 32-CPU full wrapper；测试编号 307509/307510 不是实际任务。尚未运行 MethSCAn prepare 或后续全链路。

The MethSCAn 1.1.0 CLI was verified. Its native ALLC input is used directly. The current source contains 6,554 indexed ALLCs (~19 GB; 3,165 CYL and 3,389 ZCP). A balanced four-cell intake smoke passed, validating 400 CGN records and canonical cell links. Python, shell, CLI, skill-package, and Slurm test-only checks passed; test IDs 307509/307510 are not jobs. MethSCAn prepare and the downstream chain have not yet run.

## 初版固定参数 / Preliminary fixed parameters

| Parameter | Value |
|---|---:|
| MethSCAn min sites | 60,000 |
| Mean methylation range | 20–85% |
| Smooth bandwidth | 1,000 bp |
| Scan bandwidth / step | 2,000 / 100 bp |
| Variable-window fraction | 0.02 |
| Minimum cells per VMR | 6 |
| Scan/matrix threads | 32 |
| VMR observed-cell fraction | 0.05 |
| Minimum observed VMRs/cell | 100 |
| PCA / neighbours / Leiden | 30 / 15 / 0.8 |
| Random seed | 20260825 |

这些值来自参考 `Methscan.md` 和当前 Scanpy 风格，只是首轮候选；必须由 smoke/full 结果评估后再确认。 / These values are first-pass candidates derived from `Methscan.md` and the current Scanpy conventions; smoke/full evidence is required before confirmation.

## 未决事项 / Open items

- QC 主表尚未完成，当前不指定细胞 QC 列；不得把 `NA` 当作通过。
- 参考基因组尚未确认，TSS BED/profile 默认关闭。
- 正式分析前需要计算节点读取 ALLC 源的 preflight、空间评估和 20-cell smoke。
- 还没有 VMR 数量、过滤保留率、缺失率、PCA/UMAP/Leiden 或生物学注释结果。

The QC table and genome build/TSS BED remain unresolved. A compute-node input preflight, storage estimate, and 20-cell smoke are required before a formal run. No VMR, retention, missingness, embedding, clustering, or biological annotation result exists yet.
