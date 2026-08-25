# Environment / 环境

本阶段发现、选择并验证项目所需的软件环境，支持共享环境和分阶段独立环境。

This stage discovers, selects, and validates software environments, supporting both a shared environment and isolated per-stage environments.

## 策略 / Strategy

- 优先复用已有且验证通过的环境。 / Prefer existing verified environments.
- 仅在依赖兼容时共用环境；存在版本冲突时按阶段隔离。 / Share only when dependencies are compatible; isolate stages when versions conflict.
- 未经使用者确认不创建、安装、升级或删除环境。 / Do not create, install, upgrade, or remove environments without user approval.
- 调度脚本使用已记录的绝对可执行路径。 / Scheduler scripts use recorded absolute executable paths.

最终映射维护在 `Supplementary/environments.tsv`；探测结果写入 `Results/Environment/`，日志写入 `Scripts/Environment/logs/`。

Maintain the final mapping in `Supplementary/environments.tsv`; write discovery results to `Results/Environment/` and logs to `Scripts/Environment/logs/`.

当前项目采用混合策略：兼容步骤可共享环境，存在依赖冲突的步骤使用独立环境。 / The current project uses a hybrid strategy: compatible stages share an environment, while conflicting stages remain isolated.

## 执行入口 / Entry point

```bash
bash Scripts/Environment/01_discover_environments.sh
```
