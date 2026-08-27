# MethylVI report / MethylVI 报告

## Current workflow contract / 当前流程约定

The project has two independent MethylVI routes: ALLCools-selected 5-kb bins
and MethSCAn-discovered VMRs. Both aggregate integer CGN methylated and total
coverage counts; neither trains MethylVI on transformed score matrices.

项目有两条彼此独立的 MethylVI 路线：ALLCools 筛选的 5-kb bins 与
MethSCAn 发现的 VMR。两者都使用整数 CGN `mc/cov` 计数；均不将变换后的
score 矩阵作为 MethylVI 输入。

## 2026-08-26 maintenance / 2026-08-26 维护

- VMR defaults now use `Results/Methscan/CYL_ZCP_full_20260826_final`'s
  selected-ALLC manifest and expected VMR branch.
- The upstream run has not completed its scan stage; no VMR-MethylVI run has
  been started from it.
- Slurm outputs are routed to `Scripts/Methylvi/logs/`.

- VMR 默认读取 `Results/Methscan/CYL_ZCP_full_20260826_final` 的已筛选
  ALLC manifest 和预期 VMR 分支。
- 上游运行尚未完成 scan 阶段；尚未基于它启动 VMR-MethylVI。
- Slurm 标准输出与错误输出已统一写入 `Scripts/Methylvi/logs/`。
