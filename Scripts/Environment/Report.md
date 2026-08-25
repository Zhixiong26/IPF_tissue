# Environment 分析报告 / Analysis report

状态 / Status：环境发现与初始映射已完成 / discovery and initial mapping complete

## 已验证环境 / Verified environments

- BAM/ALLCools/Scanpy：`allcools`，Samtools 1.22、ALLCools 1.1.1、Scanpy 1.9.3。
- Methscan：`MethSCAn`，MethSCAn 1.1.0、Python 3.8.20。
- MethylVI：`methylvi`，Scanpy 1.12.3、scvi-tools 1.5.0.post1、Torch 2.13.0。
- FASTQ/Bismark 候选：`moabs`，Bismark 0.24.2、Bowtie2 2.5.4；其中 Samtools 0.1.19 过旧，不用于 BAM 验证。

The selected hybrid mapping is recorded in `Supplementary/environments.tsv`. Shared environments outside this project are read-only dependencies and must not be modified.

## 发现结果 / Discovery results

- 完整只读清单位于 `Results/Environment/environment_inventory.tsv`。 / The complete read-only inventory is at `Results/Environment/environment_inventory.tsv`.
- 初次发现脚本因旧版 Bash 对空数组的 `set -u` 行为失败，已修正并重新完整运行。 / The first discovery attempt failed because of old-Bash empty-array behavior under `set -u`; it was fixed and rerun successfully.
- 未创建、安装、升级或删除任何环境。 / No environment was created, installed into, upgraded, or removed.

## 服务器资源快照 / Server resource snapshot

检查时间 / Checked: 2026-08-24 14:23 CST

- 登录节点 `ctl01`：32 个物理 CPU cores（2 × Xeon Silver 4314），125 GiB RAM，其中约 104 GiB available；swap 31 GiB，已使用约 15 GiB。登录节点仅用于轻量检查，不作为完整 QC 运行资源。 / Login node `ctl01`: 32 physical CPU cores (2 × Xeon Silver 4314), 125 GiB RAM with about 104 GiB available; 31 GiB swap with about 15 GiB used. The login node is for lightweight checks, not full QC runs.
- 项目文件系统 `/home/lijia/luozhixiong`：204 TiB 中仅约 1 GiB available，使用率 100%；inode 使用率 96%。在清理或确定新结果位置前，不应向这里写入大型中间文件。 / Project filesystem `/home/lijia/luozhixiong`: only about 1 GiB available out of 204 TiB, 100% used; inode usage is 96%. Do not write large intermediates here until space is cleared or a new output location is selected.
- 源数据盘 `/mnt/data04`：100 TiB，总体使用率 96%，约 4.7 TiB available；inode 余量充足。 / Source-data filesystem `/mnt/data04`: 100 TiB total, 96% used, about 4.7 TiB available; inode capacity is ample.
- Slurm `cpu`：3 × 56-core nodes、每节点约 257.4 GiB RAM；快照时 23/168 CPUs allocated。`cu02` 与 `cu03` 为 idle。 / Slurm `cpu`: three 56-core nodes with about 257.4 GiB RAM each; 23/168 CPUs were allocated at the snapshot. `cu02` and `cu03` were idle.
- Slurm `fat`：1 × 192-core node、约 1.03 TiB RAM；快照时 5/192 CPUs allocated，约 913 GiB free memory。 / Slurm `fat`: one 192-core node with about 1.03 TiB RAM; 5/192 CPUs were allocated and about 913 GiB memory was free at the snapshot.
- 未发现 `nvidia-smi`，当前分区也未报告 GPU GRES；本阶段按 CPU-only 规划。用户当时无 Slurm 作业。 / `nvidia-smi` was not found and current partitions reported no GPU GRES; plan this stage as CPU-only. The user had no Slurm jobs at the time.

后续复查 / Later recheck (2026-08-24)：`free -h` 中的 16 GiB 是 **free swap**，不是可用 RAM；登录节点 RAM available 约 108 GiB。项目文件系统可用空间已从约 1 GiB 回升到约 40 GiB，inode 使用率约 35%，但 `df` 仍按整数显示 100% used。QC 表格输出很小，可以运行；大型中间文件仍不应写入项目盘。

The 16 GiB shown by `free -h` is **free swap**, not available RAM; login-node available RAM is about 108 GiB. Project-filesystem free space later increased from about 1 GiB to about 40 GiB and inode usage is about 35%, although integer-rounded `df` still displays 100% used. Small QC tables can be written safely; large intermediates should still not target the project filesystem.
