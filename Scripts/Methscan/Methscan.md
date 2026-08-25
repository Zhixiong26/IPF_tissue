# MethSCAn workflow reference / MethSCAn 流程参考

Do not convert ALLC files to coverage files for this project. MethSCAn natively reads ALLC through `methscan prepare --input-format allc`; the maintained automated entry point is `run_methscan.sbatch`.

本项目不要将 ALLC 转换为 coverage 文件。MethSCAn 原生通过 `methscan prepare --input-format allc` 读取 ALLC；维护中的自动化入口为 `run_methscan.sbatch`。

## Required sequence / 固定顺序

```text
ALLC intake
→ exact Scanpy cell_id whitelist
→ exclude Scanpy cell_type = NA
→ MethSCAn prepare
→ filter: covered CpG ≥ 300,000 and overall mCG ≥ 0.50
→ smooth → VMR scan → matrix → VMR-Scanpy
```

The Scanpy whitelist is defined by `Results/Scanpy/E_CYL_ZCP_notebook/cell_id_cell_type.tsv`. The audit outputs are:

- `00_scanpy_selected/allc_excluded_by_scanpy.tsv`
- `00_scanpy_selected/input_manifest.tsv`

Scanpy `cell_type=NA` is an explicit exclusion, not a missing-value label to retain. The pre-prepare non-`NA` whitelist is the only RNA selection step. MethSCAn filter only removes cells by methylation QC, and downstream VMR-Scanpy inherits labels from the entry manifest without a second cell-type match.

Scanpy 的白名单由 `Results/Scanpy/E_CYL_ZCP_notebook/cell_id_cell_type.tsv` 的 `cell_id` 列定义。`cell_type=NA` 是明确排除条件；此筛选仅在 prepare 前执行一次，filter 后不会再次按 cell type 匹配或复核。

See `README.md` for commands and `Report.md` for verified runs and limitations.
