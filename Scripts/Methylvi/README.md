# IPF 30wcov：ALLCools → MethylVI 分析说明

本目录实现以下可复现流程：

```text
6列 Bismark coverage
  → 转换为 BGZF + tabix ALLC
  → ALLCools 5-kb CGN hypo-score MCDS
  → 二值化、区域过滤、LSI、Leiden、tSNE、UMAP、ConsensusClustering
  → 按 ALLCools 保留的细胞和 5-kb bins 从 ALLC 重建整数 mc/cov
  → MethylVI 批次校正
  → 20维 latent、15近邻、UMAP、Leiden
```

上游参数与 `git@github.com:Zhixiong26/MethylVI.git` 的 `yuanpei/` 流程一致。当前项目在其后增加 MethylVI；不会把经过处理的 ALLCools hypo-score 当成 MethylVI 原始计数。

## 服务器路径

| 项目 | 默认路径 |
|---|---|
| 仓库 | `/home/lijia/luozhixiong/IPF_tissue` |
| 脚本 | `/home/lijia/luozhixiong/IPF_tissue/Scripts` |
| Bismark coverage | `Data/30wcov` |
| coverage 实际目录 | `/home/lijia/jiangyuanpei/methscan/xunyin/IPF_tissue/30wcov` |
| 细胞类型注释 | `Supplementary/manual_celltype_annotation.tsv` |
| canonical GRCh38 sizes | `Supplementary/hg38.canonical.chrom.sizes` |
| 新流程根输出 | `Results/MethylVI_30wcov_allcools_blacklist_f0p2` |
| ALLCools 输出 | `Results/MethylVI_30wcov_allcools_blacklist_f0p2/allcools_5kb` |
| MethylVI H5MU | `Results/MethylVI_30wcov_allcools_blacklist_f0p2/ipf_allcools_5kb_methylvi_input.h5mu` |
| MethylVI 结果 | `Results/MethylVI_30wcov_allcools_blacklist_f0p2/results` |

新流程使用新的输出根目录，不会复用旧方差 top-N 流程位于 `Results/MethylVI_30wcov` 的 checkpoint。

`Data/30wcov` 和注释表是指向共享数据的符号链接。不要把约 144 GB coverage、转换后的 ALLC、MCDS、H5MU 或模型提交到 Git。

## 输入核验

2026-08-20 实测：

| 项目 | 数量 |
|---|---:|
| `*_allc.gz.cov` | 6,554 |
| 唯一 coverage cell IDs | 6,554 |
| 注释表记录 | 6,808 |
| coverage 与注释 ID 匹配 | 5,203 |
| 匹配且具有非空 `manual_celltype` | 5,111 |
| 匹配但 cell type 为空 | 92 |
| coverage 未匹配注释表 | 1,351 |

为与 `yuanpei/` 保持一致，默认将全部 6,554 个 coverage 细胞送入无监督 ALLCools 和 MethylVI。没有有效标签的细胞保留，并标记为 `manual_celltype=Unknown`。如果只想分析注释表中出现的 5,203 个细胞，可在运行前设置：

```bash
export IPF_INCLUDE_UNANNOTATED=0
export IPF_EXPECTED_CELLS=5203
```

## 两个 Conda 环境

ALLCools 与当前 scvi-tools 的 Python 版本要求不同，因此统一入口直接调用两个独立环境：

| 用途 | 环境 | 已验证版本 |
|---|---|---|
| ALLCools 上游 | `/home/lijia/jiangyuanpei/miniforge3/envs/allcools` | Python 3.9、ALLCools 1.1.1 |
| MethylVI | `/home/lijia/luozhixiong/miniconda3/envs/methylvi` | Python 3.12、scvi-tools 1.5.0.post1、PyTorch 2.13.0 CPU/MKL |

环境核验已集成到：

```bash
cd /home/lijia/luozhixiong/IPF_tissue/Scripts
bash 03_run_methylvi.sh verify
```

若共享 ALLCools 环境发生迁移，可在提交前覆盖：

```bash
export IPF_ALLCOOLS_ENV=/path/to/allcools/environment
```

## 脚本和执行阶段

| 文件 | 功能 |
|---|---|
| `00_methylvi_config.sh` | 路径、环境、ALLCools 和 MethylVI 参数 |
| `01_prepare_allcools.py` | 选择细胞，优先复用已验证 ALLC，缺失时将 `.cov` 转为 BGZF ALLC，并建立输入 manifest 和 ALLC table |
| `02_cluster_allcools.py` | `yuanpei` 风格的 hypo-score、LSI 和 consensus clustering |
| `03_build_methylvi_from_allcools.py` | 对 H5AD 保留 bins 重建整数 `mc/cov` H5MU |
| `04_train_methylvi.py` | 训练 MethylVI 并生成 latent、UMAP 和 Leiden |
| `03_run_methylvi.sh` | 统一入口和断点控制 |

统一入口：

```bash
cd /home/lijia/luozhixiong/IPF_tissue/Scripts

bash 03_run_methylvi.sh verify    # 只读核验输入、环境和细胞数
bash 03_run_methylvi.sh prepare   # cov→ALLC，并生成5-kb hypo-score MCDS
bash 03_run_methylvi.sh cluster   # 已有MCDS→ALLCools聚类H5AD
bash 03_run_methylvi.sh allcools  # prepare + cluster
bash 03_run_methylvi.sh build     # ALLCools H5AD + ALLC→整数mc/cov H5MU
bash 03_run_methylvi.sh train     # 训练MethylVI
bash 03_run_methylvi.sh smoke     # 100个平衡抽取细胞、2 epochs的独立测试
bash 03_run_methylvi.sh all       # 完整正式流程
```

`prepare`、MethylVI count-row 构建以及 MCDS/H5AD 都支持安全复用。输入细胞、保留 features 或关键参数发生变化时，manifest/config 检查会拒绝混用旧 checkpoint；应改用新的 `IPF_MVI_ROOT`，不要直接覆盖来源不一致的结果。

## 与 yuanpei 对齐的 ALLCools 参数

| 环节 | 参数 | 默认值 |
|---|---|---:|
| ALLC context | `IPF_MC_CONTEXT` | `CGN` |
| genomic bin | `IPF_BIN_SIZE` | 5,000 bp |
| MCDS quantifier | hypo-score cutoff | 0.9 |
| 二值化 | `binarize_matrix` cutoff | 0.95 |
| feature 筛选 | `filter_regions(hypo_percent=)` | 3.06%（约 50,077 bins） |
| LSI | algorithm / seed | `arpack` / 0 |
| 显著 PC | p cutoff | 0.1 |
| 图邻居 | neighbors | 25 |
| 初步 Leiden | resolution | 1.0 |
| tSNE | perplexity | 30 |
| ConsensusClustering | min cluster size | 10 |
| ConsensusClustering | Leiden repeats | 500 |
| ConsensusClustering | Leiden resolution | 0.5 |
| ConsensusClustering | consensus rate | 0.5 |
| ConsensusClustering | train fraction / max n | 0.5 / 500 |

本流程不再执行“单细胞 bin coverage ≥20、至少50细胞、方差排序前50,000”的旧特征选择。聚类前会使用 ENCODE `ENCFF356LFX` GRCh38 blacklist，移除重叠比例 ≥0.2 的 5-kb bin；最终 MethylVI features 由 blacklist 过滤和 ALLCools `filter_regions()` 共同决定。

## 计数语义

coverage 到 ALLC 的转换为：

```text
ALLC column 5 mc  = coverage column 5
ALLC column 6 cov = coverage column 5 + column 6
context            = CGN
strand             = +
```

默认优先直接引用已有的 6,554 个可读 ALLC（由同一批 `IPF_tissue/30wcov` 生成）：

```text
/home/lijia/jiangyuanpei/methscan/xunyin/IPF_tissue/allcools_5kbin/input_allc
```

文件必须同时具备非空 `.allc.tsv.gz` 和 `.tbi`。流程将路径、大小和修改时间写入
`input_allc_manifest.tsv`；只有缺失的细胞才回退到原始 `.cov` 转换。

ALLCools H5AD 的 `X` 是 hypo-score，只用于区域选择、LSI 和聚类。MethylVI builder 重新遍历逐细胞 ALLC，在保留的 5-kb bins 内汇总整数 `mc/cov`，并验证：

- cell ID 与 ALLC table 一一匹配；
- feature 坐标唯一且与 5-kb bins 对齐；
- `mc`、`cov` 为整数且始终满足 `mc <= cov`；
- checkpoint manifest 的 cell/feature SHA256 与当前 H5AD 一致；
- 根据最大 coverage 自动选择 `uint16` 或 `uint32`，绝不截断计数。

## MethylVI 参数与批次变量

| 参数 | 默认值 |
|---|---:|
| `batch_key` | `cohort`（CYL/ZCP，临时值） |
| latent | 20 |
| hidden / layers | 128 / 1 |
| likelihood | `betabinomial` |
| dispersion | `region` |
| max epochs | 500 |
| batch size | 32 |
| early stopping | 开启 |
| MethylVI neighbors | 15 |
| MethylVI Leiden resolution | 1.0 |
| seed | 0 |

重要：目前只能从 cell ID 前缀得到 `CYL/ZCP`。如果它表示疾病、组织或其他目标生物分组，就不应作为 batch，否则可能消除目标生物差异。正式解释结果前必须补充 donor/sample/technical batch 元数据，并将 `IPF_BATCH_KEY` 改成真实技术或个体批次列。

## Smoke test

```bash
cd /home/lijia/luozhixiong/IPF_tissue/Scripts
export IPF_THREADS=8
bash 03_run_methylvi.sh smoke
```

集群上应从仓库根目录提交已提供的 smoke 作业，避免在登录节点运行：

```bash
cd /home/lijia/luozhixiong/IPF_tissue
sbatch run_methylvi_allcools_smoke.sbatch
```

smoke 输出使用独立目录 `Results/MethylVI_30wcov_allcools/smoke`。默认从 CYL/ZCP 轮流抽取 100 个细胞并训练 2 epochs，因此能够验证两个 cohort 均进入流程。ALLCools 即使对少量细胞也需要建立整套 canonical GRCh38 5-kb 区间，耗时会明显长于旧的直接 `.cov` smoke test。

## SLURM 节点资源

2026-08-20 快照：

| Partition | Node | CPUs | Total memory (MiB) | Observed free memory (MiB) | State |
|---|---|---:|---:|---:|---|
| `cpu` | `cu01` | 56 | 257,400 | 27,628 | `mix` |
| `cpu` | `cu02` | 56 | 257,400 | 1,428 | `mix` |
| `cpu` | `cu03` | 56 | 257,400 | 720 | `mix` |
| `fat` | `fat01` | 192 | 1,031,600 | 938,507 | `mix` |

实时检查：

```bash
sinfo -N -o '%P|%N|%c|%m|%e|%t|%G'
squeue -p cpu,fat
```

完整流程包含 6,554 个细胞、ALLC/MCDS 构建和大规模 H5MU 组装。当前正式文件使用 `fat` 分区、`lijia` account、128 CPUs、250 GiB；该节点能够容纳正式任务的内存请求。不要在登录/控制节点运行。

2026-08-20 最初 `luozhixiong` 没有 Slurm account 关联；管理员随后已加入实验室 account：

```text
luozhixiong → account lijia → QOS normal
```

当前已通过 `sbatch --test-only` 和实际 smoke 作业验证，`lijia` account 可正常提交到 `cpu`/`fat` 分区。

## SLURM 正式提交

仓库根目录已提供 `run_methylvi_allcools.sbatch`，内容为：

```bash
#!/usr/bin/env bash
#SBATCH --job-name=ipf_allc_mvi
#SBATCH --account=lijia
#SBATCH --partition=fat
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=250G
#SBATCH --time=7-00:00:00
#SBATCH --output=/home/lijia/luozhixiong/IPF_tissue/allc_mvi_%j.out
#SBATCH --error=/home/lijia/luozhixiong/IPF_tissue/allc_mvi_%j.err

set -euo pipefail
export IPF_THREADS="$SLURM_CPUS_PER_TASK"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd /home/lijia/luozhixiong/IPF_tissue
bash Scripts/03_run_methylvi.sh all
```

提交和监控：

```bash
cd /home/lijia/luozhixiong/IPF_tissue
job_id=$(sbatch --parsable run_methylvi_allcools.sbatch)
echo "$job_id"
squeue -j "$job_id"
tail -f "allc_mvi_${job_id}.out"
```

任务结束后：

```bash
sacct -j "$job_id" --format=JobID,JobName,Partition,State,ExitCode,Elapsed,AllocCPUS,MaxRSS
```

## 主要输出

| 输出 | 内容 |
|---|---|
| `allcools_5kb/input_allc/*.allc.tsv.gz{,.tbi}` | 转换后的逐细胞 ALLC 与索引 |
| `allcools_5kb/selected_cells.allc.tsv` | ALLCools cell-to-file table |
| `allcools_5kb/mcg_5kb.mcds` | 5-kb CGN hypo-score MCDS |
| `allcools_5kb/mcg_5kb.clustered.h5ad` | 过滤 bins、LSI、Leiden、tSNE、UMAP、L1 |
| `allcools_5kb/cell_clusters.csv.gz` | ALLCools 细胞聚类表 |
| `allcools_5kb/allcools_5kb_tsne_L1.png` | ALLCools tSNE（L1） |
| `allcools_5kb/allcools_5kb_umap_L1.png` | ALLCools UMAP（L1） |
| `allcools_5kb/allcools_original_embedding_cell_type.png` | ALLCools UMAP（cell type） |
| `allcools_5kb/allcools_original_embedding_cohort.png` | ALLCools UMAP（cohort） |
| `count_rows/*.npz` | 对保留 features 聚合的逐细胞整数计数 checkpoint |
| `ipf_allcools_5kb_methylvi_input.h5mu` | MethylVI `mCG.layers['mc'/'cov']` 输入 |
| `results/model` | MethylVI 模型 |
| `results/methylvi_latent.npy` | 20维潜空间 |
| `results/methylvi_embedding.h5ad` | MethylVI neighbors、UMAP、Leiden |
| `results/methylvi_umap_cell_type.png` | MethylVI UMAP（cell type） |
| `results/methylvi_umap_condition.png` | MethylVI UMAP（CYL/ZCP batch） |
| `results/methylvi_umap_L1.png` | MethylVI UMAP（ALLCools L1） |
| `results/methylvi_umap_methylVI_leiden.png` | MethylVI UMAP（MethylVI Leiden） |

## VMR-MethylVI 独立流程

除固定 5-kb bin 分析外，项目还提供一条以 VMR（variably methylated regions）作为特征的独立 MethylVI 流程。两条流程复用相同的 6,554 个逐细胞 ALLC，但拥有不同的脚本、checkpoint 和结果目录，互不覆盖。

### 为什么不能直接使用旧 Hamming 脚本

`Data/cov/30_40_VMRs.sh` 将单 CpG 甲基化率按 0.5 二值化并计算细胞间 Hamming distance。该结果不能直接输入 MethylVI，因为 MethylVI 需要每个细胞、每个区域的整数甲基化计数 `mc` 和总覆盖度 `cov`。VMR 新流程不使用 Hamming distance，而是逐个扫描 ALLC 并在每个 VMR 内累加原始计数。

### VMR 输入审计与筛选

| 步骤 | VMR 数量 |
|---|---:|
| `Data/cov/VMR_1%.txt` 原始区间 | 37,798 |
| 去除非 canonical GRCh38 区间 | 37,539 |
| 去除 ENCODE blacklist（overlap ≥0.2） | 36,025 |

最终 BED 已按 canonical chromosome 顺序和起始坐标排序，区间互不重叠，并添加 `VMR_000001` 形式的唯一 ID：

```text
Results/MethylVI_30wcov_vmrs_blacklist_f0p2/input/vmr_canonical_blacklist_f0p2.bed
```

### VMR 计数逻辑

对每个细胞执行以下操作：

1. 顺序读取 BGZF ALLC，只保留 `CGN` context。
2. 将 ALLC 的 1-based CpG 位置转换成 BED 使用的 0-based 坐标。
3. 判断 CpG 所属的非重叠 VMR。
4. 在 VMR 内累加 ALLC 第5列作为 `mc`、第6列作为 `cov`。
5. 将每个细胞的稀疏计数保存为独立 `.npz` checkpoint。

默认仅保留在超过 200 个细胞（约 3.06%）中满足 `cov > 0` 的 VMR。该阈值仅影响最终矩阵组装，不影响逐细胞 checkpoint；修改阈值后可以复用已有计数。

最终输入为：

```text
Results/MethylVI_30wcov_vmrs_blacklist_f0p2/ipf_vmr_methylvi_input.h5mu
```

其中：

```text
mCG.layers["mc"]   # 细胞 × VMR 的整数甲基化计数
mCG.layers["cov"]  # 细胞 × VMR 的整数总 coverage
```

同时保存 `manual_celltype`、`cohort`、VMR 坐标、长度以及每个 VMR 的 covered-cell 数。根据最大 coverage 自动选择 `uint16` 或 `uint32`，避免计数溢出。

### VMR-MethylVI 参数

| 参数 | 默认值 |
|---|---:|
| batch key | `cohort`（CYL/ZCP） |
| likelihood | `betabinomial` |
| dispersion | `region` |
| latent dimensions | 20 |
| hidden / layers | 128 / 1 |
| epochs | 最大 500，early stopping |
| batch size | 32 |
| neighbors | 15 |
| Leiden resolution | 1.0 |

训练后生成普通 UMAP，以及 target weight 为 `0.2、0.5、0.7、0.9` 的细胞类型监督 UMAP。`Unknown` 标签以 `-1` 传入 UMAP，不作为已知细胞类型进行监督。

### VMR 脚本和运行命令

独立说明和脚本位于：

```text
VMR_MethylVI/README.md
VMR_MethylVI/Scripts/
```

分阶段运行：

```bash
cd /home/lijia/luozhixiong/IPF_tissue
bash VMR_MethylVI/Scripts/04_run_vmr_methylvi.sh verify
bash VMR_MethylVI/Scripts/04_run_vmr_methylvi.sh prepare
bash VMR_MethylVI/Scripts/04_run_vmr_methylvi.sh build
bash VMR_MethylVI/Scripts/04_run_vmr_methylvi.sh train
bash VMR_MethylVI/Scripts/04_run_vmr_methylvi.sh plots
bash VMR_MethylVI/Scripts/04_run_vmr_methylvi.sh supervised
```

完整 Slurm 作业：

```bash
sbatch VMR_MethylVI/run_methylvi_vmrs.sbatch
```

当前配置为 `cpu` 分区、56 CPUs、120 GiB 内存；首次正式运行作业为 `307457`。

主要输出：

| 输出 | 内容 |
|---|---|
| `input/prepare_summary.json` | canonical/blacklist/ALLC 输入审计 |
| `count_rows/*.npz` | 逐细胞 VMR mc/cov checkpoint |
| `build_summary.json` | coverage 分布、最终 VMR 数、dtype 和最大 coverage |
| `ipf_vmr_methylvi_input.h5mu` | VMR-MethylVI 整数计数输入 |
| `results/model` | MethylVI 模型 |
| `results/methylvi_latent.npy` | 20维 MethylVI latent |
| `results/methylvi_embedding.h5ad` | neighbors、UMAP、Leiden 和元数据 |
| `results/methylvi_vmr_umap_cell_type.png` | 按细胞类型着色的 UMAP |
| `results/methylvi_vmr_umap_condition.png` | 按 CYL/ZCP batch 着色的 UMAP |
| `results/methylvi_vmr_umap_methylVI_leiden.png` | 按 MethylVI Leiden 着色的 UMAP |

## 修改记录

| 日期 | 状态 | 修改 | 验证 |
|---|---|---|---|
| 2026-08-20 | 本地待提交 | 将直接 `.cov` 方差 top-N 流程改为 `yuanpei` 风格 ALLCools 上游，再重建整数计数训练 MethylVI | 两个环境导入、6,554 文件审计、Bash/Python 语法、真实 `.cov→ALLC→tabix`、2细胞 MCDS、2×100 H5MU 构建均通过 |
| 2026-08-22 | 运行中 | 新增 canonical + blacklist VMR-MethylVI 独立流程；36,025 个候选 VMR，按 covered-cell 数进一步过滤 | Bash/Python 语法、6,554 ALLC 输入、2个真实细胞 VMR mc/cov H5MU smoke test 均通过；正式作业 `307457` |

正式任务完成后，应在本表追加 Git commit、SLURM job ID、最终细胞数、保留 bin 数、训练停止 epoch 和 `MaxRSS`。
