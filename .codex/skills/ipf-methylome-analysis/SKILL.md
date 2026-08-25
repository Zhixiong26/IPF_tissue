---
name: ipf-methylome-analysis
description: Iteratively analyse this IPF tissue project across methylome QC, transcriptome Scanpy/Harmony, ALLCools, methscan, and MethylVI. Use when planning, running, checking, documenting, or extending these project workflows.
---

# IPF tissue multi-omic analysis

Work with the user as an iterative analyst: inspect the current project state, make the smallest justified analysis change, verify its outputs, and record observed facts before proceeding. Do not claim that a stage is complete from a submitted job alone.

The project root is `/home/lijia/luozhixiong/IPF_tissue`. Read [the project map](references/project-map.md) before choosing a stage and [the analysis state](references/analysis-state.md) at the start of each session. Update the analysis state after a material observation, completed run, or change to the workflow. Keep it factual: include commands/job IDs, inputs, outputs, versions, QC summaries, and unresolved decisions.

## Required project layout

At the beginning of a new analysis session, run `scripts/ensure_project_layout.sh` or equivalently ensure these project-root directories exist:

- `Scripts/`: versioned analysis scripts and scheduler entry points.
- `Results/`: stage- and run-specific analysis products; never place raw input here.
- `Supplementary/`: sample metadata, annotations, reference tables, and other small supporting inputs.

Keep the capitalization exactly as shown. The legacy root-level `log/` and `logs/` directories may contain historical outputs; preserve them, but do not route new logs there.

Organize new code and outputs into matching stage directories under both `Scripts/` and `Results/`. Each `Scripts/<stage>/` must also contain its own `logs/` directory; all stdout, stderr, scheduler logs, and run diagnostics for that stage go there. The initial stages are `Environment/`, `FastQ/`, `QC/`, `Scanpy/`, `Methscan/`, and `Methylvi/`. Add another matching pair when the workflow gains a substantial new stage; pass additional stage names after the project path to `scripts/ensure_project_layout.sh`. Existing root-level scripts and historical result directories remain in place unless the user explicitly requests a migration.

Treat `Data/` as the local input/staging area and `Results/` as the derived-output area. Do not infer compute-node accessibility from a project symlink: the current fat node cannot access `/mnt/data04`, so resolve and validate the actual target from the intended execution context before submitting. Never route a compute job directly to `/mnt/data04` for this project. Stage required inputs beneath the project only when the user has authorized the copy and storage has been checked.

Every `Scripts/<stage>/` directory must contain both `README.md` and `Report.md`. For this project, maintain both documents bilingually in Chinese and English until the user changes this convention. Keep the two languages semantically aligned in the same sections. The README is the maintained operational guide: purpose, inputs, outputs, environments, script order, and entry points. The report is the evidence ledger: executed commands/jobs, parameters, QC, verified results, failures, decisions, and limitations. Update both as the stage evolves. `scripts/ensure_project_layout.sh` creates missing bilingual documents but preserves existing files.

## GitHub collaboration gate

At the start of a new project, remind the user to create a GitHub repository and provide its SSH clone URL, expected default branch, and whether the local project should be initialized or connected to that repository. If this information is missing, prompt:

> 请先建立并提供 GitHub 仓库的 SSH 地址（例如 `git@github.com:OWNER/REPO.git`），并说明默认分支。请同时确认当前服务器是否已配置 GitHub SSH；不要发送私钥、密码或访问令牌。

Inspect existing local Git state and remotes before changing them. Creating a remote repository, changing a remote, committing, or pushing requires the user's explicit request. If SSH is not configured, guide the user through local key generation and adding only the **public key** to GitHub; never print, request, store, or commit a private key or token. Before the first commit, ensure raw FASTQs, large intermediate/results files, model checkpoints, and runtime logs are excluded from Git tracking. Read [GitHub and SSH setup](references/github-ssh.md) when repository setup or troubleshooting is needed.

## Analysis starting-point gate

After the GitHub collaboration gate, ask whether the analysis should start from Raw FASTQ or skip FASTQ and begin from an existing downstream dataset. Do not assume FASTQ is required, and do not treat files left over from another run as authorization to use them.

If the starting point is not already explicit, pause and prompt:

> 本项目是否需要从 Raw_fastq 开始？如果需要，请提供 Raw_fastq 路径，并说明是否要软链接到 `Data/Raw_fastq/`；如果跳过 FASTQ，请提供实际起始数据类型（例如 BAM、coverage、ALLC 或矩阵）、路径和已知的生成流程/参考基因组信息。

If the user chooses FASTQ, verify that the path exists and is readable, then ask only for still-missing facts such as paired-end versus single-end, library protocol, and sample/group metadata. Create links only when requested or clearly agreed; never copy or move large raw data by default. Record the confirmed source path and link mapping in `references/analysis-state.md`.

If the user skips FASTQ, record the supplied downstream source and known provenance, defer FASTQ work without deleting any existing artifacts, and validate the selected starting point before use. For BAM, confirm readability, file/index pairing, BAM quick-check, headers, reference sequence dictionary, read groups, sort order, and sample/group mapping. Apply format-appropriate checks for coverage, ALLC, matrices, or other inputs. Do not infer that skipped upstream QC or processing passed merely because downstream files exist.

## Environment discovery and configuration

After selecting the starting data, discover existing executables and Conda/Mamba environments before proposing installation. Use `Scripts/Environment/01_discover_environments.sh` and store its inventory under `Results/Environment/`. Ask whether the user prefers one shared environment or separate stage environments when that preference is not already explicit.

Support both modes:

- **Shared:** reuse one environment only when all selected stages have compatible dependencies and a smoke import/command check passes.
- **Per-stage:** assign independent environments to FastQ/BAM-Bismark, Methscan, Scanpy/ALLCools, MethylVI, or other stages when requirements conflict or isolation improves reproducibility.

Prefer existing verified environments. Do not create, install into, upgrade, or remove an environment without the user's approval. Never silently modify a shared environment owned by another user. Record chosen absolute environment and executable paths in the relevant bilingual README/Report and maintain the project mapping in `Supplementary/environments.tsv`. Scheduler scripts should call absolute executables or use an explicit recorded activation method; do not rely on an interactive shell's active environment. Verify imports and tool versions with a small smoke check before full analysis.

## Stage routing

1. **Project and input intake.** Apply the GitHub collaboration gate, then the analysis starting-point gate. For FASTQ, use `scripts/inspect_fastq_layout.sh`. When FASTQ is skipped, validate the user-supplied downstream dataset first. Derive explicit sample metadata from evidence; do not infer read structure, cell barcodes, or biological conditions from directory names alone.
2. **Environment.** Discover existing tools/environments, choose shared versus per-stage assignment based on compatibility and user preference, record absolute paths, and smoke-test before installing or running analyses.
3. **Bismark.** Before alignment, establish the library protocol, genome build and Bismark genome directory. Perform/read FastQC or equivalent QC, then use protocol-appropriate trimming, alignment, methylation extraction, and coverage generation. Keep per-sample logs and mapping/conversion/duplication metrics. Do not deduplicate reads unless that matches the library protocol.
4. **Per-cell QC.** Keep FASTQ input pairs, BAM primary mapped pairs/records, and ALLC methylation metrics as separately recoverable layers, then join on explicit sample ID and cell barcode. For the linked Met FASTQs, the verified structure is R1 bases 1–17 as the cell barcode and bases 18–29 as the UMI; the 8-bp FASTQ-header suffix is the sample index. Generate the trusted barcode set from BAM `CB` tags first and use the step-02 TSV as the step-01 whitelist. Count BAM pairs using primary mapped read1 records carrying `CB`; do not use BAM-internal mapped percentage as mapping efficiency. Use bounded file-level process parallelism for the scan-heavy steps: FASTQ by partition, BAM by prefix file, and ALLC by cell; keep small TSV joins serial. Allocate resources by actual parallel units: the current Met data have two FASTQ partitions but 156 prefix BAMs per sample, so use the separate `07_run_parallel_bam.sbatch` and dependent `08_run_parallel_fastq.sbatch` arrays rather than reserving BAM-sized CPU allocations during FASTQ scans. Treat absent ALLC contexts as missing, never as zero methylation. Preserve all per-metric flags and do not delete failed cells during the first pass. Use `Scripts/QC/01_FASTQ_per_cell_reads.py` through `06_QC_summary.py` and document project-specific thresholds separately from primary QC.
5. **methscan.** Confirm the locally installed methscan command/version and the requested biological operation before adding an executable invocation. Preserve its exact input/output contract in the state file; do not substitute a similarly named tool or invent flags.
6. **Transcriptome Scanpy.** For CYL/ZCP transcriptome-E QC, Scrublet, Harmony, clustering, annotation, plotting, export, or submission, read [scanpy-transcriptome.md](references/scanpy-transcriptome.md). For cluster validation, label correction, or marker-dotplot changes, also read [scanpy-cluster-validation.md](references/scanpy-cluster-validation.md). Use the copied matrices in `Data/Matrix/`, the validated notebook, and the production scripts under `Scripts/Scanpy/Scripts/`. Keep CYL/ZCP cohort metadata separate from cell type, exclude LC from this analysis, and never reuse a cluster-number mapping when the observed cluster set changes.
7. **Methylome Scanpy / ALLCools.** For the existing single-cell coverage workflow, use `Scripts/03_run_methylvi.sh` and the configuration it sources. Inspect the resulting H5AD metadata, dimensions, QC, neighbours and embeddings before interpreting clusters. Do not confuse this coverage-derived workflow with the CYL/ZCP transcriptome-E Scanpy workflow.
8. **MethylVI.** Use the integer methylated-count and coverage layers built from ALLC/coverage, not the ALLCools hypo-score matrix. Validate that counts are integral and `mc <= cov`, retain input manifests, and treat checkpoint reuse as valid only when its cell/feature/config signatures match.

## Operational rules

- Put new analysis outputs below the matching `Results/<stage>/` directory, using a distinct named run root, and put execution logs under `Scripts/<stage>/logs/`. Never overwrite an existing result root, raw FASTQs, shared coverage, or a model checkpoint merely to rerun a stage.
- Prefer a small, representative smoke test before a full scheduler submission. Show the exact command and requested compute resources before a costly full run; wait for the user's go-ahead when a choice of reference, protocol, or compute scale is still unresolved.
- Use read-only inspection to diagnose failures. Change scripts/configuration only when the evidence identifies the needed change, then validate the smallest affected stage.
- Produce tabular machine-readable QC/manifests alongside figures. Record sample identity and condition in explicit metadata rather than relying on directory names downstream.
- Keep every production script's explicit and fixed analysis parameters synchronized with its stage `Report.md`, including QC, normalization, PCA/integration, clustering, marker/annotation, plotting, output schema, environment, and scheduler resources. Record runtime-derived values after execution rather than presenting them as pre-run facts.
- Existing `Scripts/` and `VMR_MethylVI/` files are prior work. Preserve them unless the current task specifically requires a compatible edit.

## References

- Read [project-map.md](references/project-map.md) for local paths, existing entry points, and known environment boundaries.
- Read and update [analysis-state.md](references/analysis-state.md) throughout this iterative project.
- Read [scanpy-transcriptome.md](references/scanpy-transcriptome.md) for CYL/ZCP RNA/Scanpy work.
- Read [scanpy-cluster-validation.md](references/scanpy-cluster-validation.md) for cluster evidence review, annotation correction, and merged cell-type dotplots.
- Read [github-ssh.md](references/github-ssh.md) only when establishing or troubleshooting GitHub/SSH collaboration.
